immutable Cosmological <: SpaceType
    a3::Float64
    H2::Float64
    adda::Float64
end

function Cosmological(t=0.001)
    const a = A(t)
    Cosmological(a*a*a, 2*H(t), ADDA(t))
end
