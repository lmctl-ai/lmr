# The discrete Fourier transform: fft(1:8) and the amplitude spectrum of a two-tone signal.
source("r/_emit.R")
z <- fft(1:8)
t <- 0:63
s <- sin(2 * pi * 5 * t / 64) + 0.5 * cos(2 * pi * 12 * t / 64)
emit("fft-spectrum",
     fft_1_8_re = Re(z), fft_1_8_im = Im(z), spectrum_mod = Mod(fft(s))[1:33],
     inverse_roundtrip = Re(fft(fft(s), inverse = TRUE) / 64))
