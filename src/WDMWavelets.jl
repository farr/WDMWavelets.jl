module WDMWavelets

using FFTW, SpecialFunctions

export Phi_unit, fd_wavelet_basis_matrix, td_wavelet_basis_matrix, wdm_transform, wdm_dT_dF, wtm_times_frequencies

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

@doc raw"""
    C_matrix([T,] nt, nf)

Construct the ``C`` matrix from Eq. (10) of [Cornish
(2020)](https://link.aps.org/doi/10.1103/PhysRevD.102.124038).

The ``C`` matrix is used to pick out alternately the imaginary or real parts of
the Fourier transform to construct the wavelet decomposition.  It is alternately
`1` or `1im` according to whether the sum of the indices is even or odd.

The ``C`` matrix is defined in the same way for the zero and highest frequency
components (`C[:,1]` and `C[:,nf]`), but the extraction of the wavelet
decomposition should not use it for these special cases; see [Necula, Klimenko,
& Mitselmakher
(2012)](https://iopscience.iop.org/article/10.1088/1742-6596/363/1/012032) for
details.
"""
function C_matrix(T, nt, nf)
    C = zeros(T, nt, nf)
    for i in axes(C, 1)
        for j in axes(C, 2)
            if (i+j-2) % 2 == 0
                C[i, j] = 1
            else
                C[i, j] = 1im
            end
        end
    end
    C
end

"""
    fd_wavelet_basis_matrix(nt, nf, A, d)

Returns the fourier domain representation of the wavelet basis functions.

Returns a complex tensor of size `(nt, nf, nt*nf)` that represents the Fourier
domain representation of the wavelet basis elements.
"""
function fd_wavelet_basis_matrix(nt, nf, A, d)
    n = nt*nf

    dT, dF = wdm_dT_dF(nt, nf, 1)

    C = C_matrix(nt, nf)

    ns = 0:nt-1
    ms = 0:nf-1
    fs = fftfreq(n)

    ns = reshape(ns, nt, 1, 1)
    ms = reshape(ms, 1, nf, 1)
    fs = reshape(fs, 1, 1, n)
    C = reshape(C, nt, nf, 1)

    @. exp(-1im * ns * 2 * pi * fs * dT) * (C * Phi_unit(fs / dF - ms, A, d) + conj(C) * Phi_unit(fs / dF + ms, A, d)) / sqrt(2*dF)
end

"""
    td_wavelet_basis_matrix(nt, nf, A, d)

Returns the time-domain representation of the wavelet basis functions.

Returns a real tensor of size `(nt, nf, nt*nf)` that represents the time domain
representation of the wavelet basis elements; literally the inverse Fourier
transform of the Fourier domain representation.

These basis functions are orthonormal (see the associated test).

TODO: currently the lowest-frequency and highest-frequency basis functions are
incorrectly computed presently, and are not orthonormal.
"""
function td_wavelet_basis_matrix(nt, nf, A, d)
    real(ifft(fd_wavelet_basis_matrix(nt, nf, A, d), 3))
end

"""
    wdm_dT_dF(nt, nf, dt)

Returns the time and frequency bin widths of the WDM transform.
"""
wdm_dT_dF(nt, nf, dt) = (nf*dt, 1/(2*nf*dt))

"""
    wdm_times_frequencies(nt, nf, dt)

Returns `(ts, fs)` giving the times and frequencies of the corresponding columns
and rows of the wdm matrix.  `nt` and `nf` are the number of time and frequency
bins, and `dt` is the sample rate of the input signal.
"""
function wdm_times_frequencies(nt, nf, dt)
    dT = nf*dt
    dF = 1/(2*dT)

    (dT * (0:nt-1), dF * (0:nf-1))
end

"""
    wdm_transform(x, nt, nf, A, d)

Returns the matrix of WDM coefficients for the input signal `x`.

The matrix will have `nt` time bins and `nf` frequency bins; the wavelet
parameters are `A` and `d` (see `Phi_unit` above).  

The returned matrix will be of size `(nt, nf)`.
"""
function wdm_transform(x, nt, nf, A, d)
    n = nt*nf
    @assert nt % 2 == 0
    @assert nf % 2 == 0
    @assert length(x) == n

    nto2 = div(nt, 2)

    _, dF = wdm_dT_dF(nt, nf, 1)

    X = fft(x)
    X = fftshift(X)

    fs = fftfreq(n)
    fs_phi = vcat(fs[1:nto2], fs[end-nto2+1:end])
    phi = @. Phi_unit(fs_phi / dF, A, d) / sqrt(dF)
    phi = fftshift(phi)
    xmn = zeros(Complex{Float64}, nt, nf)
    for i in axes(xmn, 2)
        i0 = div(n,2) + 1 + (i-1)*nto2
        xmn[:,i] = ifft( ifftshift(X[i0-nto2:i0+nto2-1] .* phi) ) * nt
    end

    C = C_matrix(nt, nf)
    sign_matrix = [ (i*j % 2 == 0 ? 1 : -1) for i in 0:nt-1, j in 0:nf-1 ]
    result = @. sqrt(2) * sign_matrix * real(C * xmn) / n
    result
end

function wdm_inverse_transform(x, A, d)
    nt, nf = size(x)

    n = nt * nf

    nto2 = div(nt, 2)
    no2 = div(n,2)

    _, dF = wdm_dT_dF(nt, nf, 1)
    
    fs = fftfreq(n)
    fs_phi = vcat(fs[1:nto2], fs[end-nto2+1:end])
    phi = @. Phi_unit(fs_phi / dF, A, d) / sqrt(dF)
    phi = fftshift(phi)

    xf = zeros(Complex{Float64}, n)

    xc = zeros(Complex{Float64}, nt, nf)
    for i in axes(xc, 1)
        for j in axes(xc, 2)
            if (i + j - 2) % 2 == 0
                xc[i,j] = x[i,j]
            elseif (j-1) % 2 == 0
                xc[i,j] = -1im * x[i,j]
            else
                xc[i,j] = 1im * x[i,j]
            end
        end
    end
    xc = fft(xc, 1)

    for j in axes(xc, 2)
        i0 = div(n,2) + 1 + (j-1)*nto2
        xf[i0-nto2:i0+nto2-1] += fftshift(xc[:,j]) .* phi
    end
    xf = ifftshift(xf)
    for i in 0:no2-1
        xf[end-i] = conj(xf[i+2])
    end

    real.(ifft(xf)) ./ sqrt(2)
end 

end # module WDMWavelets
