# TODO: parametrize with real cosmological quantities like Ω_m, etc.

# for a matter dominated universe:

# as function of real time
@inline A(t) = 0.0225*t^(2/3)

# as function of expansion coefficient `a`
@inline Ta(a) = 296.296296*a*sqrt(a)

###################

@inline FD(a1,a2) = 20.0*(1./sqrt(a1)-1./sqrt(a2))
@inline FK(a1,a2) = 20.0*(sqrt(a2) - sqrt(a1))

@inline Ha(a) = 0.1/a/a*sqrt(a)
