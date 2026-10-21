// Tests for rcli::fft: against a naive O(n^2) DFT in long double for many lengths, and against R 4.3.3 fft() output.
#include "rcli/fft.hpp"

#include <cstdint>
#include <string>

#include "check.hpp"

using cd = std::complex<double>;
using rcli::CVec;

namespace {
std::uint64_t state = 88172645463325252ULL;
double unif() {  // xorshift64, uniform in [-1, 1)
  state ^= state << 13;
  state ^= state >> 7;
  state ^= state << 17;
  return static_cast<double>(state >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}

// Naive DFT in long double at the output bins listed in ks.
CVec naive(const CVec& z, int sign, const std::vector<size_t>& ks) {
  const size_t n = z.size();
  const long double pi = 3.141592653589793238462643383279502884L;
  std::vector<long double> c(n), s(n);  // root table so the O(n^2) loop has no trig calls
  for (size_t j = 0; j < n; ++j) {
    const long double a = 2 * pi * static_cast<long double>(j) / static_cast<long double>(n);
    c[j] = std::cos(a);
    s[j] = sign * std::sin(a);
  }
  CVec out(ks.size());
  for (size_t ik = 0; ik < ks.size(); ++ik) {
    const size_t k = ks[ik];
    long double re = 0, im = 0;
    size_t e = 0;  // (j*k) mod n
    for (size_t j = 0; j < n; ++j) {
      re += z[j].real() * c[e] - z[j].imag() * s[e];
      im += z[j].real() * s[e] + z[j].imag() * c[e];
      e += k;
      if (e >= n) e -= n;
    }
    out[ik] = cd(static_cast<double>(re), static_cast<double>(im));
  }
  return out;
}

void check_vec(const CVec& got, const std::vector<cd>& want, double tol, const std::string& what) {
  check::ok(got.size() == want.size(), what + " length");
  for (size_t i = 0; i < want.size() && i < got.size(); ++i) {
    check::close(got[i].real(), want[i].real(), tol, what + " Re[" + std::to_string(i + 1) + "]");
    check::close(got[i].imag(), want[i].imag(), tol, what + " Im[" + std::to_string(i + 1) + "]");
  }
}
}  // namespace

int main() {
  // Naive DFT comparison: error <= 1e-12 * max|z| * log2(n+1), both directions.
  double worst = 0;
  for (size_t n : {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 16, 25, 97, 100, 101, 128, 211, 1000, 1009, 4096, 30030, 2 * 1009}) {
    CVec z(n);
    for (auto& v : z) v = cd(unif(), unif());
    for (int sign : {-1, 1}) {
      // all bins for n <= 1024, else 1..n/64-strided bins plus a few at the ends (the naive DFT is O(n) per bin)
      std::vector<size_t> ks;
      const size_t stride = n <= 1024 ? 1 : n / 64;
      for (size_t k = 0; k < n; k += stride) ks.push_back(k);
      if (stride > 1) for (size_t k : {size_t(1), size_t(2), n / 2, n - 1}) ks.push_back(k);
      const CVec got = rcli::fft(z, sign > 0);
      const CVec want = naive(z, sign, ks);
      double maxz = 0, err = 0;
      for (const auto& v : z) maxz = std::fmax(maxz, std::abs(v));
      for (size_t i = 0; i < ks.size(); ++i) err = std::fmax(err, std::abs(got[ks[i]] - want[i]));
      const double bound = 1e-12 * maxz * std::log2(static_cast<double>(n) + 1.0);
      worst = std::fmax(worst, err / (maxz * std::log2(static_cast<double>(n) + 1.0)));
      check::ok(err <= bound, "fft n=" + std::to_string(n) + (sign > 0 ? " inverse" : " forward") +
                                  " err " + std::to_string(err) + " bound " + std::to_string(bound));
    }
  }
  std::fprintf(stderr, "fft: worst err / (max|z| log2(n+1)) = %.3g\n", worst);
  // Round trip: fft(fft(z), inverse=TRUE)/n == z for a Bluestein length
  {
    CVec z(1009);
    for (auto& v : z) v = cd(unif(), unif());
    const CVec back = rcli::fft(rcli::fft(z), true);
    double err = 0;
    for (size_t i = 0; i < z.size(); ++i) err = std::fmax(err, std::abs(back[i] / 1009.0 - z[i]));
    check::ok(err < 1e-13, "fft round trip n=1009");
  }

  // z <- complex(real=1:7, imaginary=c(0,-1,0.5,2,0,0,1)); fft(z)
  const CVec z7 = {{1, 0}, {2, -1}, {3, 0.5}, {4, 2}, {5, 0}, {6, 0}, {7, 1}};
  check_vec(rcli::fft(z7),
            {{28, 2.5}, {-3.7084315306100293, 5.3546266852201834}, {-7.2304606588584841, 3.5876520308546711},
             {-2.8088273951054825, 0.66555519338226421}, {-4.1911726048945139, -0.93214912734878774},
             {0.23046065885848765, -1.9946616913221558}, {-3.2915684693899672, -9.1810230907861747}},
            1e-14, "fft(z7)");
  // fft(z, inverse=TRUE)
  check_vec(rcli::fft(z7, true),
            {{28, 2.5}, {-3.2915684693899672, -9.1810230907861747}, {0.23046065885848765, -1.9946616913221558},
             {-4.1911726048945139, -0.93214912734878774}, {-2.8088273951054825, 0.66555519338226421},
             {-7.2304606588584841, 3.5876520308546711}, {-3.7084315306100293, 5.3546266852201834}},
            1e-14, "fft(z7, inverse=TRUE)");
  // x <- cos((0:100)/3) + (0:100)/50; fft(x)[c(1,2,3,51,101)]   (n = 101 prime -> Bluestein)
  {
    rcli::Vec x(101);
    for (int i = 0; i <= 100; ++i) x[i] = std::cos(i / 3.0) + i / 50.0;
    const CVec f = rcli::fft_real(x);
    const CVec pick = {f[0], f[1], f[2], f[50], f[100]};
    check_vec(pick,
              {{104.12555849261146, 0}, {2.1997398934820382, 33.413635798586192},
               {2.4931666323924384, 18.350498590213927}, {-0.26104902769912514, 0.0026841393183088453},
               {2.1997398934820382, -33.413635798586192}},
              1e-13, "fft(cos((0:100)/3) + (0:100)/50)");
  }
  // x <- sin((0:11)^2/7); fft(x)
  {
    rcli::Vec x(12);
    for (int i = 0; i < 12; ++i) x[i] = std::sin(i * i / 7.0);
    check_vec(rcli::fft_real(x),
              {{1.15797105614338, 0}, {0.20622459071072519, -1.8576210835681679},
               {-2.6199359857542177, -1.1180015719358081}, {0.41229733307467037, 1.7297266943814829},
               {-1.7607744942821744, -2.7203196600075992}, {2.1077895366286241, -1.8077501918574095},
               {2.1508269831013664, 0}, {2.1077895366286228, 1.8077501918574099},
               {-1.760774494282175, 2.7203196600075992}, {0.41229733307467159, -1.7297266943814829},
               {-2.6199359857542182, 1.1180015719358081}, {0.20622459071072408, 1.8576210835681675}},
              1e-14, "fft(sin((0:11)^2/7))");
  }
  check::ok(rcli::fft(CVec{}).empty(), "fft(complex(0)) is empty");
  check::close(rcli::fft(CVec{{3, -2}})[0].imag(), -2, 0, "fft(3-2i) == 3-2i");
  check::throws([] { rcli::fft_real({1, rcli::kNaN}); }, "fft with NaN errors");
  return check::report("test_fft");
}
