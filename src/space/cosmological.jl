
immutable Cosmological <: SpaceType
    F1::Float64   # === 1 / (  sa * a^5  )
    F2::Float64   # === ddt(dadt) / a^3 / sa
    F3::Float64   # === 3./2a
    sa::Float64
end

function Cosmological(a=0.01)
    a3 = a*a*a
    sa = Sa(a)

    Cosmological(
        1./sa/a3/a/a,
        ADDa(a)/a3/sa,
        3./2a,
        sa
    )
end
