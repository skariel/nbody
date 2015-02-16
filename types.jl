immutable Particle <: AbstractPoint3D
    _x::Float64
    _y::Float64
    _z::Float64
    _m::Float64
end
Particle() = Particle(0., 0., 0., 0.)
getx(p::Particle) = p._x
gety(p::Particle) = p._y
getz(p::Particle) = p._z

type World
    tree::OctTree{Particle}
    particles::Array{Particle, 1}
    vx::Array{Float64, 1}
    vy::Array{Float64, 1}
    vz::Array{Float64, 1}
    ax::Array{Float64, 1}
    ay::Array{Float64, 1}
    az::Array{Float64, 1}
    n::Int64
    opening_alpha2::Float64
    opening_excluded_frac2::Float64
    smth2::Float64
    dtfrac::Float64
end
