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
    tree::CompiledOctTree{Particle}
    particles::SharedArray{Particle, 1}
    vx::SharedArray{Float64, 1}
    vy::SharedArray{Float64, 1}
    vz::SharedArray{Float64, 1}
    ax::SharedArray{Float64, 1}
    ay::SharedArray{Float64, 1}
    az::SharedArray{Float64, 1}
    n::Int64
    opening_alpha2::Float64
    smth2::Float64
    dtfrac::Float64
end

type Simulation
    tree::OctTree{Particle}
    w::World
    t::Float64
    dt::Float64
    ti::Float64
    tf::Float64
    stepf::Int64
    xi::SharedArray{Float64, 1}
    yi::SharedArray{Float64, 1}
    zi::SharedArray{Float64, 1}
    vxi::SharedArray{Float64, 1}
    vyi::SharedArray{Float64, 1}
    vzi::SharedArray{Float64, 1}
    step::Int64
    limit_by_steps::Bool
end

function Simulation(w::World; ti=0.0, tf=1.0, stepc=100, limit_by_steps=false)
    xi = SharedArray(Float64, w.n)
    yi = SharedArray(Float64, w.n)
    zx = SharedArray(Float64, w.n)
    vxi = SharedArray(Float64, w.n)
    vyi = SharedArray(Float64, w.n)
    vzi = SharedArray(Float64, w.n)
    for i in 1:w.n
        @inbounds xi[i] = w.particles[i]._x
        @inbounds yi[i] = w.particles[i]._y
        @inbounds zi[i] = w.particles[i]._z
        @inbounds vxi[i] = 0.0
        @inbounds vyi[i] = 0.0
        @inbounds vzi[i] = 0.0
    end
    Simulation(
        createtree(w),
        w,
        ti,  # t
        0.0, # dt
        ti,  # ti
        tf,  # tf
        stepc,   # stepf::Int64
        xi, yi, zi, vxi, vyi, vzi,
        0, # step::Int64
        limit_by_steps # limit_by_steps::Bool
    )
end

function reset!(s::Simulation)
    # set initial positions and velocities
    for i in 1:s.w.n
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
    tree::Array{CompiledOctTree{Particle}, 1}
    dt::Array{Float64, 1}
    w::World
    steps::Int64
end

type Optimization
    history::History
    x0::SharedArray{FLoat64,1}
    y0::SharedArray{FLoat64,1}
    z0::SharedArray{FLoat64,1}
    gx::SharedArray{FLoat64,1}
    gy::SharedArray{FLoat64,1}
    gz::SharedArray{FLoat64,1}
    function Optimization(sim::Simulation)
        hist = History(CompiledOctTree{Particle}[], Float64[], sim.w, 0)
        x0 = SharedArray(Float64, sim.w.n)
        y0 = SharedArray(Float64, sim.w.n)
        z0 = SharedArray(Float64, sim.w.n)
        gx = SharedArray(Float64, sim.w.n)
        gy = SharedArray(Float64, sim.w.n)
        gz = SharedArray(Float64, sim.w.n)
        for i in 1:sim.w.n
            @inbounds x0[i] = sim.xi[i]
            @inbounds y0[i] = sim.yi[i]
            @inbounds z0[i] = sim.zi[i]
            @inbounds gx[i] = 0.0
            @inbounds gy[i] = 0.0
            @inbounds gz[i] = 0.0
        end
        new(hist, x0,y0,z0, gx,gy,gz)
    end
end
