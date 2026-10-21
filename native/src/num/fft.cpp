// Discrete Fourier transform of arbitrary length (clean-room).
// Method:
//   * Mixed-radix decimation-in-time Cooley-Tukey (Cooley & Tukey 1965; Singleton 1969 "An algorithm for computing
//     the mixed radix fast Fourier transform", IEEE Trans. Audio Electroacoust. 17, 93-103) with butterflies for
//     radix 4, 2, 3, 5 and a direct O(p^2) DFT butterfly for any other prime factor p.
//   * Bluestein's chirp-z algorithm (Bluestein 1970; Rabiner, Schafer & Rader 1969) when n has a prime factor > 100:
//     X_k = conj(w_k) * sum_j (z_j conj(w_j)) w_{k-j}, w_m = exp(i*pi*m^2/n), evaluated as a circular convolution with
//     a power-of-two FFT.
// Twiddle factors are computed directly from long double sin/cos of an exactly reduced angle (index mod n), so the
// error grows like O(eps log n).
#include "rcli/fft.hpp"

#include <cstdint>

namespace rcli {
namespace {

using cd = std::complex<double>;
constexpr long double kPiL = 3.141592653589793238462643383279502884L;

// exp(sign * 2*pi*i * k / n) for 0 <= k < n, computed in long double with octant symmetry left to libm.
cd unit_root(std::uint64_t k, std::uint64_t n, int sign) {
  const long double a = 2.0L * kPiL * static_cast<long double>(k) / static_cast<long double>(n);
  return cd(static_cast<double>(std::cos(a)), static_cast<double>(sign * std::sin(a)));
}

class Plan {
 public:
  Plan(size_t n, int sign) : n_(n), sign_(sign), tw_(n) {
    for (size_t k = 0; k < n; ++k) tw_[k] = unit_root(k, n, sign);
    size_t m = n;
    while (m % 4 == 0) { factors_.push_back(4); m /= 4; }
    while (m % 2 == 0) { factors_.push_back(2); m /= 2; }
    for (size_t p = 3; p * p <= m; p += 2)
      while (m % p == 0) { factors_.push_back(p); m /= p; }
    if (m > 1) factors_.push_back(m);
  }

  void run(const cd* in, cd* out) const { rec(in, 1, out, n_, 0); }

 private:
  size_t n_;
  int sign_;
  CVec tw_;
  std::vector<size_t> factors_;

  void rec(const cd* in, size_t stride, cd* out, size_t n, size_t level) const {
    if (n == 1) { out[0] = in[0]; return; }
    const size_t p = factors_[level];
    const size_t m = n / p;
    for (size_t q = 0; q < p; ++q) rec(in + q * stride, stride * p, out + q * m, m, level + 1);
    const size_t step = n_ / n;  // W_n^e = W_N^{e*step}
    const size_t pstep = n_ / p;  // W_p^e = W_N^{e*pstep}
    std::vector<cd> t(p), u(p);
    for (size_t k = 0; k < m; ++k) {
      t[0] = out[k];
      for (size_t q = 1; q < p; ++q) t[q] = out[q * m + k] * tw_[q * k * step];
      switch (p) {
        case 2:
          u[0] = t[0] + t[1];
          u[1] = t[0] - t[1];
          break;
        case 4: {
          const cd a = t[0] + t[2], b = t[0] - t[2], c = t[1] + t[3];
          const cd d0 = t[1] - t[3];
          const cd d = sign_ < 0 ? cd(d0.imag(), -d0.real()) : cd(-d0.imag(), d0.real());  // sign*i*(t1-t3)
          u[0] = a + c;
          u[2] = a - c;
          u[1] = b + d;
          u[3] = b - d;
          break;
        }
        case 3: {
          const double s3 = sign_ * 0.86602540378443864676;  // sign*sin(2pi/3)
          const cd sum = t[1] + t[2], dif = t[1] - t[2];
          const cd mid = t[0] - 0.5 * sum;
          const cd rot(-s3 * dif.imag(), s3 * dif.real());  // i*s3*dif
          u[0] = t[0] + sum;
          u[1] = mid + rot;
          u[2] = mid - rot;
          break;
        }
        case 5: {
          const double c1 = 0.30901699437494742410, c2 = -0.80901699437494742410;  // cos(2pi/5), cos(4pi/5)
          const double s1 = sign_ * 0.95105651629515357212, s2 = sign_ * 0.58778525229247312917;
          const cd a1 = t[1] + t[4], b1 = t[1] - t[4], a2 = t[2] + t[3], b2 = t[2] - t[3];
          const cd r1 = t[0] + c1 * a1 + c2 * a2, r2 = t[0] + c2 * a1 + c1 * a2;
          const cd i1 = s1 * b1 + s2 * b2, i2 = s2 * b1 - s1 * b2;
          const cd ii1(-i1.imag(), i1.real()), ii2(-i2.imag(), i2.real());
          u[0] = t[0] + a1 + a2;
          u[1] = r1 + ii1;
          u[4] = r1 - ii1;
          u[2] = r2 + ii2;
          u[3] = r2 - ii2;
          break;
        }
        default:
          for (size_t s = 0; s < p; ++s) {
            cd acc = t[0];
            for (size_t q = 1; q < p; ++q) acc += t[q] * tw_[((q * s) % p) * pstep];
            u[s] = acc;
          }
      }
      for (size_t s = 0; s < p; ++s) out[s * m + k] = u[s];
    }
  }
};

bool has_large_prime(size_t n) {
  size_t m = n;
  for (size_t p = 2; p * p <= m; ++p)
    while (m % p == 0) {
      if (p > 100) return true;
      m /= p;
    }
  return m > 100;
}

CVec bluestein(const CVec& z, int sign) {
  const size_t n = z.size();
  size_t m = 1;
  while (m < 2 * n - 1) m <<= 1;
  // chirp c_j = exp(sign*i*pi*j^2/n) = W_{2n}^{j^2 mod 2n} (exact integer reduction of the angle)
  CVec chirp(n);
  const std::uint64_t two_n = 2 * static_cast<std::uint64_t>(n);
  for (size_t j = 0; j < n; ++j) {
    const std::uint64_t jj = static_cast<std::uint64_t>(j);
    chirp[j] = unit_root((jj * jj) % two_n, two_n, sign);
  }
  CVec a(m, cd(0.0, 0.0)), b(m, cd(0.0, 0.0));
  for (size_t j = 0; j < n; ++j) a[j] = z[j] * chirp[j];
  // X_k = c_k * sum_j a_j conj(c_{k-j}), since jk = (j^2 + k^2 - (k-j)^2)/2
  b[0] = std::conj(chirp[0]);
  for (size_t j = 1; j < n; ++j) b[j] = b[m - j] = std::conj(chirp[j]);
  Plan fwd(m, -1), inv(m, +1);
  CVec fa(m), fb(m), conv(m);
  fwd.run(a.data(), fa.data());
  fwd.run(b.data(), fb.data());
  for (size_t k = 0; k < m; ++k) fa[k] *= fb[k];
  inv.run(fa.data(), conv.data());
  CVec out(n);
  const double scale = 1.0 / static_cast<double>(m);
  for (size_t k = 0; k < n; ++k) out[k] = chirp[k] * conv[k] * scale;
  return out;
}

}  // namespace

CVec fft(const CVec& z, bool inverse) {
  for (const cd& v : z)
    if (std::isnan(v.real()) || std::isnan(v.imag())) throw Error("fft: missing values (NaN) are not allowed");
  const size_t n = z.size();
  if (n <= 1) return z;
  const int sign = inverse ? +1 : -1;
  if (has_large_prime(n)) return bluestein(z, sign);
  Plan plan(n, sign);
  CVec out(n);
  plan.run(z.data(), out.data());
  return out;
}

CVec fft_real(const Vec& x, bool inverse) {
  CVec z(x.size());
  for (size_t i = 0; i < x.size(); ++i) z[i] = cd(x[i], 0.0);
  return fft(z, inverse);
}

}  // namespace rcli
