module WDMWavelets

using FFTW, SpecialFunctions

export Phi_unit, wdm_transform, wdm_inverse_transform, wdm_dT_dF, wdm_times_frequencies

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

    # We are going to ignore the first frequency bin because it's fucked up!
    # For the other frequency bins, we have:
    # wnm = sum_k x[k] * gnm[k]
    # wnm = sum_k (sum_l exp(2*pi*i*l*k/N) X[l]) (sum_j exp(2*pi*i*j*k/N) Gnm[j])
    # wnm = sum_l X[l] Gnm[N-l]
    # wnm = sum_l X[l] exp(-2*pi*i*n*(N-l)*Nf/N) (Cnm*G[(N-l) - m*Nt/2] + conj(Cnm)*G[(N-l) + m*Nt/2]))
    # wnm = sum_l X[l] exp(2*pi*i*n*l*Nf/N) (Cnm*G[l + m*Nt/2] + conj(Cnm)*G[l - m*Nt/2])
    # wnm = sum_l X[l] exp(2*pi*i*n*l/Nt) (Cnm*G[l + m*Nt/2] + conj(Cnm)*G[l - m*Nt/2])
    # wnm = 2*Re conj(Cnm) sum_l exp(2*pi*i*n*l/Nt) X[l] G[l - m*Nt/2]

    fs = fftfreq(n)
    fs_phi = vcat(fs[1:nto2], fs[end-nto2+1:end])
    phi = @. Phi_unit(fs_phi / dF, A, d) / sqrt(dF)

    xnm = zeros(Complex{Float64}, nt, nf)

    # Load up the columns of xmn:
    for m in 2:nf
        l0 = (m-1)*div(nt,2) + 1
        xnm[1:div(nt,2),m] .= X[l0:l0+div(nt,2)-1] .* phi[1:div(nt,2)] # Positive frequencies
        xnm[div(nt,2)+1:end,m] .= X[l0-div(nt,2):l0-1] .* phi[div(nt,2)+1:end] # Negative frequencies
    end 
    ifft!(xnm, 1)  # In-place inverse FFT

    result = zeros(nt, nf)
    for n in 1:nt
        for m in 2:nf
            C = ((n+m-2) % 2 == 0) ? 1 : 1im
            result[n,m] = sqrt(2) * real(conj(C) * xnm[n,m]) / nf # Normalization is * (nt / n), which is *nf.
        end
    end
    result
end

# x[k] = sum_{nm} wnm gnm[k]
# x[k] = sum_{nm} wnm sum_l exp(2*pi*i*l*k/N) Gnm[l]
# x[k] = sum_{nm} wnm sum_l exp(2*pi*i*l*k/N) exp(-2*pi*i*n*l/Nt) (Cnm * G[l - m*Nt/2] + conj(Cnm) * G[l + m*Nt/2])])
# x[k] = sum_l sum_m exp(2*pi*i*l*k/N) (G[l - m*Nt/2] * Ylm + G[l + m*Nt/2] * conj(Ylm))
# => X[l] = (G[l - m*Nt/2] * Ylm + G[l + m*Nt/2] * conj(Y(-l)m))
function wdm_inverse_transform(x, A, d)
    nt, nf = size(x)

    n = nt * nf

    nto2 = div(nt, 2)

    _, dF = wdm_dT_dF(nt, nf, 1)
    
    fs = fftfreq(n)
    fs_phi = vcat(fs[1:nto2], fs[end-nto2+1:end])
    phi = @. Phi_unit(fs_phi / dF, A, d) / sqrt(dF)

    ylm = zeros(Complex{Float64}, nt, nf)
    for n in 1:nt
        for m in 2:nf
            C = ((n+m-2) % 2 == 0) ? 1 : 1im
            ylm[n,m] = C * x[n,m] / sqrt(2)
        end
    end
    fft!(ylm, 1)

    X = zeros(Complex{Float64}, n)
    for m in 2:nf
        l0 = (m-1)*div(nt,2) + 1
        X[l0:l0+div(nt,2)-1] .+= ylm[1:div(nt,2),m] .* phi[1:div(nt,2)] # Positive frequencies
        X[l0-div(nt,2):l0-1] .+= ylm[div(nt,2)+1:end,m] .* phi[div(nt,2)+1:end] # Negative frequencies

        # The FFT of the conjugate is the backwards of the FFT of the original.
        l1 = n - (m-1)*div(nt,2) + 1
        X[l1] = conj(ylm[1,m]) * phi[1] # Zero frequency

        X[l1+1:l1+div(nt,2)-1] .+= conj(ylm[end:-1:div(nt,2)+2,m]) .* phi[2:div(nt,2)] # Positive frequencies
        X[l1-div(nt,2):l1-1] .+= conj(ylm[div(nt,2)+1:-1:2,m]) .* phi[div(nt,2)+1:end] # Negative frequencies
    end
    ifft!(X)
    real.(X)
end 

end # module WDMWavelets
