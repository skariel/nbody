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
withxyz(p::Particle, x::Float64, y::Float64, z::Float64) = Particle(x,y,z, p._m)
addxyz(p::Particle, dx::Float64, dy::Float64, dz::Float64) = Particle(p._x+dx,p._y+dy,p._z+dz, p._m)

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

type Simulation
    w::World
    t::Float64
    dt::Float64
    ti::Float64
    tf::Float64
    stepf::Int64
    xi::Array{Float64, 1}
    yi::Array{Float64, 1}
    zi::Array{Float64, 1}
    vxi::Array{Float64, 1}
    vyi::Array{Float64, 1}
    vzi::Array{Float64, 1}
    step::Int64
    limit_by_steps::Bool
end

Simulation(w::World; ti=0.0, tf=1.0, stepc=100, limit_by_steps=false) =
    Simulation(
        w,
        ti,  # t
        0.0, # dt
        ti,  # ti
        tf,  # tf
        stepc,   # stepf::Int64
        [p._x for p in w.particles], #xi
        [p._y for p in w.particles], #yi
        [p._z for p in w.particles], #zi
        [0.0  for p in w.particles], #vxi
        [0.0  for p in w.particles], #vyi
        [0.0  for p in w.particles], #vzi
        0, # step::Int64
        limit_by_steps # limit_by_steps::Bool
    )

function reset!(s::Simulation)
    # set initial positions and velocities
    @inbounds for i in 1:s.w.n
        s.w.particles[i] = withxyz(s.w.particles[i], s.xi[i], s.yi[i], s.zi[i])
        s.w.vx[i] = 0.0
        s.w.vy[i] = 0.0
        s.w.vz[i] = 0.0
        s.vxi[i] = 0.0
        s.vyi[i] = 0.0
        s.vzi[i] = 0.0
    end
    # set times and steps
    s.step = 0
    s.t = s.ti
end

type History
    tree::Array{OctTree{Particle}, 1}
    dt::Array{Float64, 1}
    w::World
end
