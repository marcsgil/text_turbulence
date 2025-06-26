using CairoMakie, FFTW

N = 1024
L = 1.0
Δx = L / N
xs = (0:N-1) * Δx
c = 1.0

u0 = exp.(-((xs .- 0.5 * L) ./ 0.1) .^ 2)
u̇0 = 2c * xs .* u0 / 2
ks = fftfreq(N, 2π / Δx)

u̇0 ./ ks

val = u0 .- im * u̇0 ./ ks / c
val[1] = 0.0  # Set the DC component to zero to avoid division by zero

αs = fft(val)
βs = fft(conj(val))

ts = (0:N-1) * Δx / c
result = zeros(N, N)

for n ∈ axes(result, 2)
    result[:, n] = real(bfft(circshift(αs, n - 1)) +
                        bfft(circshift(βs, 1 - n)))
end

lines(result[:, 9])