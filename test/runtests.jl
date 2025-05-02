using Test, WDMWavelets

A_wavelet = 0.25
d_wavelet = 4

function chirp(ts, Ac, As, f, fdot)
    phases = @. 2 * pi * ts * (f + fdot * ts / 2)
    @. Ac * cos(phases) + As * sin(phases)
end

@testset "WDMWavelets.jl Tests" begin
    @testset "Wavelet Transform of Chirp" begin
        dt = 1/pi
        fny = 1/(2*dt)

        nt = 64
        nf = 64
        n = nt*nf

        ts = dt * (0:n-1)
        T = n*dt

        f0 = fny / 10
        fdot = f0/T

        A = 1
        phi = atan(randn(), randn())
        Ac = A * cos(phi)
        As = A * sin(phi)

        f = chirp(ts, Ac, As, f0, fdot)

        f_tilde = wdm_transform(f, nt, nf, A_wavelet, d_wavelet)

        @testset "Parseval's Theorem" begin
            @test isapprox(sum(f.^2), sum(f_tilde.^2), rtol=1e-3, atol=0)
        end

        @testset "Chirp Track" begin
            dT, dF = wdm_dT_dF(nt, nf, dt)

            max_power_index = [x[2] for x in argmax(abs.(f_tilde), dims=2)]
            predicted_max_power_index = (f0 .+ fdot .* (0:nt-1) .* dT)/dF

            @assert all(abs.(max_power_index[2:end-1] .- predicted_max_power_index[2:end-1]) .<= 2)
        end
    end

    @testset "Testing Single Element Inverse" begin
        nt = 32
        nf = 32

        for i in 1:nt
            for j in 2:nf
                x = zeros(nt, nf)
                x[i,j] = 1.0

                xx = wdm_transform(wdm_inverse_transform(x, A_wavelet, d_wavelet), nt, nf, A_wavelet, d_wavelet)

                @test all(isapprox.(x, xx, atol=1e-8))
            end
        end
    end
end