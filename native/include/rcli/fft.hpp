// Discrete Fourier transform of any length (mixed radix 2/3/5/generic, Bluestein for large primes).
// Unnormalised in both directions, exactly like R fft(z, inverse=): forward uses exp(-2*pi*i*j*k/n),
// inverse uses exp(+2*pi*i*j*k/n) and does NOT divide by n. Accuracy: within 1e-12 * max|z| * log2(n) of a naive DFT.
#pragma once
#include <complex>
#include <vector>
#include "rcli/common.hpp"

namespace rcli {

using CVec = std::vector<std::complex<double>>;
CVec fft(const CVec& z, bool inverse = false);
CVec fft_real(const Vec& x, bool inverse = false);

}  // namespace rcli
