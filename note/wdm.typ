#set page(
  paper: "us-letter",
  numbering: "1",
)
#set par(justify: true)
#set text(
  font: "Times New Roman",
  size: 12pt,
)

#align(center, text(17pt)[
  The WDM Wavelet Basis
])

#align(center)[
    #link("https://orcid.org/0000-0003-1540-8562")[Will M. Farr] \
    #link("mailto:will.farr@stonybrook.edu") \
    #link("mailto:wfarr@flatironinstitute.org") \
    Department of Physics and Astronomy, Stony Brook University, Stony Brook NY 11794, USA \
    Center for Computational Astrophysics, Flatiron Institute, New York NY 10010, USA 
]

Some notes deriving the forward and inverse transforms used in the code.  See @Cornish2020 or @Necula2012 for full details.

The underlying "mother wavelet" for the WDM basis is compact in the frequency domain, has parameters $A$, $B$, and $d$, and is defined by 
$ tilde(phi)(f | A, B, d) #sym.prop cases(
  1 "," quad & |f| < A,
  cos(pi / 2 nu_d ((|f|-A)/B)) "," quad & A <= |f| < A + B,
  0 "," quad & |f| >= A + B
), $
where $nu_d$ is the normalized incomplete Beta function,
$ nu_d (x) = (integral_0^x  y^(d-1) (1-y)^(d-1) dif y ) / (integral_0^1 y^(d-1) (1-y)^(d-1) dif y) . $
For a bandwidth $Delta F$, $A$ and $B$ satisfy $2A + B = Delta F$; $d$ controls the sharpness of the decay to zero at frequencies between $A$ and $A+B$.  

Given a choice of the number of time bins $N_t$, and number of frequency bins $N_f$, with $N_t N_f = N$, the $n m$ ($n = 0, dots, N_t - 1$, $m = 0, dots, N_f$) component of the wavelet basis in the Fourier domain is given by 
$ tilde(g)_(n m)(f) := cases(
  e^(-pi i n f Delta T) tilde(phi)(f) "," quad & n = 0,
  e^(-2 pi i n f Delta T) (C_(n m) tilde(phi)(f - m Delta F) + C_(n m)^* tilde(phi)(f + m Delta F)) ", " quad & n > 0
) $
with 
$ C_(n m) = cases(
  1 "," quad & n + m " even",
  i "," quad & n + m " odd"
), $
and $Delta T Delta F = 1/2$.  Each basis element corresponds to a version of the mother wavelet time-shifted by $n Delta T$ ($m > 0$) or $2 n (Delta T)$ ($m = 0$) and centered at frequencies $ plus.minus m Delta F$.  For a uniformly sampled data segment of length $N$, over a time $T$, we have the definitions
$ 
delta t & = T/N \
delta f & = 1/T \
Delta T & = N_f delta t \
Delta F & = 1/(2 Delta T) = (N_t / 2) delta f,
$
so that the parameters $N_t$ and $N_f$ scale the frequency and time resolution of the wavelet basis.  

The treatment of the zero-frequency ($m = 0$) and Nyquist frequency ($m = N_f$) bands is special---and I don't quite understand it.  Currently the code does not handle the $m = 0$ wavelet band at all (enforces zeros in both directions).

The forward transform can be derived by writing 
$
  w_(n m) & = sum_k x_k g_(n m)(k delta t) \
    & = sum_k (1 / N^2) (sum_l tilde(x)_l e^(2 pi i (k l) / N)) (sum_j tilde(g)_(n m)(j delta f) e^(2 pi i (k j) / N)).
$
Performing the sum over $k$ gives $N$ times a delta function enforcing $j + l = N$; summing over $j$ yields $j = N-l$, and we have 
$
  w_(n m) = sum_l tilde(x)_l tilde(g)_(n m)((N-l) delta f).
$
Inserting the defition of $tilde(g)_(n m)$ in terms of the mother wavelet gives 
$
  w_(n m) = sum_l e^(-2 pi i n ((N-l) delta f) Delta T) tilde(x)_l (C_(n m) tilde(phi)((N-l) delta f - m Delta F) + C_(n m)^* tilde(phi)((N-l) delta f + m Delta F)).
$
Simplifying the exponential, recalling that $delta f Delta T = 1/(N_t)$ gives 
$
  w_(n m) = sum_l e^(2 pi i (n l)/N_t) tilde(x)_l (C_(n m) tilde(phi)((N-l) delta f - m Delta F) + C_(n m)^* tilde(phi)((N-l) delta f + m Delta F))
$
The arguments of the frequency domain mother wavelets simplify to
$
  (N-l) delta f - m Delta F = (N - l - m N_t / 2) delta f tilde (l + m N_t/2) delta f
$
and
$
  (N-l) delta f + m Delta F = (N - l + m N_t / 2) delta f tilde (l - m N_t/2) delta f,
$
where we have used the fact that, in a Fourier transform of the mother wavelet with frequency resolution $delta f$, the coefficients at index $i$ and $N-i$ are equal to express the final equivalence.  

Thus we have 
$
  w_(n m) = sum_l e^(2 pi i (n l)/N_t) tilde(x)_l (C_(n m) tilde(phi)_(l + m N_t/2) + C_(n m)^* tilde(phi)_(l - m N_t/2)),
$
where $tilde(phi)_k$ is the $k$th coefficient of the length-$N$ Fourier transform of the mother wavelet.  Conveniently, the exponential is periodic in $l$ with period $N_t$, and due to the band-limited structure of the mother wavelet, $tilde(phi)_k$ is non-zero only for $0 <= k <= N_t/2$ and $N - N_t/2 < k < N$.

Since the transform is real, we are free to compute the transform of one of these terms and then double the real part, so we have 
$
  w_(n m) = 2 bb(R) C_(n m)^* y_(n m),
$
where 
$
  y_(n m) = sum_l e^(2 pi i (n l)/N_t) tilde(x)_l tilde(phi)_(l - m N_t/2).
$
The quantity $y_(n m)$ can be computed for each $m$ by inverse Fourier transforming the length $N_t$ slice of $tilde(x)$ centered at a zero frequency index at $l = m N_t/2$ extending from $(m-1) N_t /2 < l <= (m+1) N_t / 2$ against the corresponding frequencies of the mother wavelet and storing the result into the $n$ entries.  There are $N_f$ such length-$N_t$ FFTs (again, ignoring the $m = 0$ band, which must be treated differently).

The inverse transform is given similarly, starting from 
$
  x_k & = sum_(n m) w_(n m) g_(n m)(k delta t) \
  & = sum_(n m) w_(n m) 1/N sum_l tilde(g)_(n m)(l delta f) e^(2 pi i (k l) / N) \
  & = sum_(n m) w_(n m) 1/N sum_l e^(2 pi i (k l) / N) e^(-2 pi i (n l)/N_t) (C_(n m) tilde(phi)_(l-m N_t/2) + C_(n m)^* tilde(phi)_(l + m N_t /2)) \
  & = sum_l 1/N e^(2 pi i (k l) / N) sum_m (tilde(phi)_(l - m N_t/2) tilde(v)_(l m) + tilde(phi)_(l + m N_t/2) tilde(v)_((N-l) m)),
$
where 
$
  tilde(v)_(l m) = sum_n e^(-2 pi i (n l)/N_t) w_(n m) C_(n m)
$
is the Fourier transform of $w$ on $C$ over the $n$ index.  Because $w_(n m)$ is real, the $C_(n m)^*$ term is the complex conjugate of the $C_(n m)$ term, and therefore 
$
  sum_n e^(-2 pi i (n l)/N_t) w_(n m) C^*_(n m) = tilde(v)_((N-l) m).
$
Let 
$
  tilde(x)_l = sum_m (tilde(phi)_(l - m N_t/2) tilde(v)_(l m) + tilde(phi)_(l + m N_t/2) tilde(v)_((N-l) m));
$
$tilde(x)$ can be computed by packing the results of $m$ length-$N_t$ Fourier transforms that produce $tilde(v)$ into the appropriate locations of the $l$ index.  Then the time-domain signal is given by 
$
  x_k = 1/N sum_l e^(2 pi i (k l) / N) tilde(x)_l,
$
a single length-$N$ inverse Fourier transform.


#bibliography("wdm.bib", style: "american-physics-society")

