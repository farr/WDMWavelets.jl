module WDMWavelets

using CairoMakie, FFTW, SpecialFunctions

@doc raw"""
    Phi_unit(f, A, d)

Meyer window function for the WDM wavelet transform.

See Eq. (10) of [Cornish
(2020)](https://link.aps.org/doi/10.1103/PhysRevD.102.124038).  The frequency
`f` and half-width of the flat part of the window, `A`, are in units of the
frequency spacing ``\Delta f``.  
"""
function Phi_unit(f, A, d)
    B = 1 - 2*A
    if abs(f) < A
        one(f)
    elseif abs(f) < A + B
        p, _ = beta_inc(d, d, (abs(f)-A)/B)
        cos(pi * p / 2)
    else
        zero(f)
    end
end

function C_matrix(nt, nf)
    C_matrix(Complex{Float64}, nt, nf)
end

function C_matrix(T, nt, nf)
    C = zeros(T, nt, nf)
    for i in axes(C, 1)
        for j in axes(C, 2)
            if (i+j) % 2 == 0
                C[i, j] = 1
            else
                C[i, j] = 1im
            end
        end
    end
    C
end

function fd_wavelet_basis_matrix(nt, nf, dt, A, d)
    n = nt*nf
    T = n*dt

    dT = nf*dt
    dF = 1/(2*dT)

    C = C_matrix(nt, nf)

    ns = 0:nt-1
    ms = 0:nf-1
    fs = fftfreq(n, 1/dt)

    ns = reshape(ns, nt, 1, 1)
    ms = reshape(ms, 1, nf, 1)
    fs = reshape(fs, 1, 1, n)
    C = reshape(C, nt, nf, 1)

    @. exp(-1im * ns * 2 * pi * fs * dT) * (C * Phi_unit(fs / dF - ms, A, d) + conj(C) * Phi_unit(fs / dF + ms, A, d)) / sqrt(2*dF)
end

function wdm_transform(x, dt, nt, nf, A, d)
    n = nt*nf
    @assert nt % 2 == 0
    @assert nf % 2 == 0
    @assert length(x) == n

    nto2 = div(nt, 2)

    T = n*dt

    dT = nf*dt
    dF = 1/(2*dT)

    X = fft(x)

    fs = fftfreq(n, 1/dt)
    fs_phi = vcat(fs[1:nto2], fs[end-nto2+1:end])
    phi = @. Phi_unit(fs_phi / dF, A, d) / sqrt(2*dF) # sqrt(dF) ensures that square integrated phi == 1.
    xmn = zeros(Complex{Float64}, nt, nf)
    for i in axes(xmn, 2)
        Xs = circshift(X, -(i-1)*nto2)
        Xsub = vcat(Xs[1:nto2], Xs[end-nto2+1:end])
        xmn[:,i] = ifft( Xsub .* phi ) * nt
    end

    C = C_matrix(nt, nf)
    sign_matrix = [ (i*j % 2 == 0 ? 1 : -1) for i in 0:nt-1, j in 0:nf-1 ]
    result = @. sqrt(2) * sign_matrix * real(C * xmn)
    result
end

function chirp(ts, Ac, As, f, fdot)
    phases = @. 2 * pi * ts * (f + fdot * ts / 2)
    @. Ac * cos(phases) + As * sin(phases)
end

end # module WDMWavelets
