# TODO: parametrize with real cosmological quantities like Ω_m, etc.

# for a matter dominated universe:

# as function of real time
@inline A(t) = t^(2/3)
@inline ADD(t) = -2/9*t^(-4/3) # === d/dt (da/dt)
@inline H(t) = 2/3/t
@inline Z(t) = 1./A(t)-1.
@inline S(t) = 2/3

# as function of expansion coefficient `a`
@inline Ta(a) = a*sqrt(a)
@inline ADDa(a) = -2/9/a/a
@inline Ha(a) = 2/3/a^(3/2)
@inline Za(a) = 1./a-1.
@inline Sa(a) = 2/3

###################

FD(a1,a2) = 3*(1./sqrt(a1)-1./sqrt(a2))
FK(a1,a2) = 3*(sqrt(a2) - sqrt(a1))
