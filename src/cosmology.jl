# TODO: parametrize with real cosmological quantities like Ω_m, etc.

# for a matter dominated universe:

# as function of real time
@inline A(t) = t^(2/3)
@inline AD(t) = 2/3*t^(-1/3)   # === da/dt
@inline ADD(t) = -2/9*t^(-4/3) # === d/dt (da/dt)
@inline H(t) = 2/3/t
@inline Z(t) = 1./A(t)-1.

# as function of expansion coefficient `a`
@inline Ta(a) = a*sqrt(a)
@inline ADa(a) = 2/3/sqrt(a)
@inline ADDa(a) = -2/9/a/a
@inline Ha(a) = 2/3/a^(3/2)
@inline Za(a) = 1./a-1.
