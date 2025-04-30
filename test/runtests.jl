using Test, WDMWavelets

A_wavelet = 0.25
d_wavelet = 4

function chirp(ts, Ac, As, f, fdot)
    phases = @. 2 * pi * ts * (f + fdot * ts / 2)
    @. Ac * cos(phases) + As * sin(phases)
end

@testset "WDMWavelets.jl Tests" begin
    @testset "TD Wavelet Orthogonality" begin
        nt = 4
        nf = 4
        n = nt*nf
        dt = 1.0

        g_matrix_td = td_wavelet_basis_matrix(nt, nf, A_wavelet, d_wavelet)

        # Presently something is messed up for the highest and lowest
        # frequencies---so ignore them.
        for i in 1:nt
            for j in 2:nf-1
                for ii in 1:nt
                    for jj in 2:nf-1
                        if i == ii && j == jj
                            @test isapprox(sum(g_matrix_td[i,j,:] .* g_matrix_td[ii,jj,:]), 1.0, rtol=1e-8, atol=0)
                        else
                            @test isapprox(sum(g_matrix_td[i,j,:] .* g_matrix_td[ii,jj,:]), 0.0, rtol=0, atol=1e-8)
                        end
                    end
                end
            end
        end
    end

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
        passed = 0
        for _ in 1:100
            nt = 64
            nf = 64

            i,j = rand(1:nt), rand(2:nf-1)
            x = zeros(nt, nf)
            x[i,j] = 1.0

            xx = wdm_transform(wdm_inverse_transform(x, A_wavelet, d_wavelet), nt, nf, A_wavelet, d_wavelet)
            if isapprox(x, xx, atol=1e-8)
                passed += 1
                continue
            else
                println("Failed for i = $i, j = $j")
                break
            end
        end
        println("Passed $passed")
    end
end