# TODO: parametrize with real cosmological quantities like Ω_m, etc.

@inline A(t) = t^(2/3)
@inline T(a) = a*sqrt(a)
@inline AD(t) = 2/3*t^(-1/3)
@inline ADD(t) = -2/9*t^(-4/3)
@inline ADDA(t) = -2/9/t/t
@inline H(t) = 2/3/t
@inline Z(t) = 1./A(t)-1.
