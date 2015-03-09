immutable Cosmological <: SpaceType
    a3::Float64      # === a^3
    ddaa3::Float64   # === d/dt(da/dt) / a^3
    FU::Float64      # === (2H da/dt + d/dt(da/dt)) / a^2
end

function Cosmological(a=0.01)
    a2 = a*a
    a3 = a2*a
    add = ADDa(a)
    Cosmological(
        a3,
        add/a3,
        (2Ha(a)*ADa(a) + add)/a2
    )
end
