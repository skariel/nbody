immutable Particle <: AbstractPoint3D
    _x::Float64
    _y::Float64
    _z::Float64
    _m::Float64
end

Particle() = Particle(0., 0., 0., 0.)
Particle(x::Float64, y::Float64, z::Float64) = Particle(x,y,z,0.0)

getx(p::Particle) = p._x
gety(p::Particle) = p._y
getz(p::Particle) = p._z

# stuff to ease handling with immutability :)
withxyz(p::Particle, x::Float64, y::Float64, z::Float64) = Particle(x,y,z, p._m)
addxyz(p::Particle, dx::Float64, dy::Float64, dz::Float64) = Particle(p._x+dx,p._y+dy,p._z+dz, p._m)
