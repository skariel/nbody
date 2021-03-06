# TODO: initial positions for test particle
# TODO: enable execution w/o test particles for performance reasons
type Simulation{T<:SpaceType}
    tree::OctTree{Particle}
    w::World{T}
    t::Float64 # time parameter: real time for newtonian space, scale factor `a` for cosmological space
    ti::Float64
    dt::Float64
    tf::Float64
    stepf::Int64
    xi::SharedArray{Float64, 1}
    yi::SharedArray{Float64, 1}
    zi::SharedArray{Float64, 1}
    vxi::SharedArray{Float64, 1}
    vyi::SharedArray{Float64, 1}
    vzi::SharedArray{Float64, 1}
    test_particle_x::SharedArray{Float64, 1}
    test_particle_y::SharedArray{Float64, 1}
    test_particle_z::SharedArray{Float64, 1}
    test_particle_vx::SharedArray{Float64, 1}
    test_particle_vy::SharedArray{Float64, 1}
    test_particle_vz::SharedArray{Float64, 1}
    test_particle_ax::SharedArray{Float64, 1}
    test_particle_ay::SharedArray{Float64, 1}
    test_particle_az::SharedArray{Float64, 1}
    step::Int64
    limit_by_steps::Bool
end

function Simulation{T<:SpaceType}(w::World{T}; ti=0.0, tf=1.0, stepc=100, limit_by_steps=false, n_test_particle=0)
    xi = SharedArray(Float64, w.n)
    yi = SharedArray(Float64, w.n)
    zi = SharedArray(Float64, w.n)
    vxi = SharedArray(Float64, w.n)
    vyi = SharedArray(Float64, w.n)
    vzi = SharedArray(Float64, w.n)
    for i in 1:w.n
         xi[i] = w.particles[i]._x
         yi[i] = w.particles[i]._y
         zi[i] = w.particles[i]._z
         # TODO: Zeldovich approx for cosmological sims
         vxi[i] = 0.0
         vyi[i] = 0.0
         vzi[i] = 0.0
    end
    test_particle_x  = SharedArray(Float64, n_test_particle)
    test_particle_y  = SharedArray(Float64, n_test_particle)
    test_particle_z  = SharedArray(Float64, n_test_particle)
    test_particle_vx = SharedArray(Float64, n_test_particle)
    test_particle_vy = SharedArray(Float64, n_test_particle)
    test_particle_vz = SharedArray(Float64, n_test_particle)
    test_particle_ax = SharedArray(Float64, n_test_particle)
    test_particle_ay = SharedArray(Float64, n_test_particle)
    test_particle_az = SharedArray(Float64, n_test_particle)

    sim = Simulation{T}(
        createtree(w),
        w,
        ti,  # t
        ti,  # ti
        0.0, # dt
        tf,  # tf
        stepc,   # stepf::Int64
        xi, yi, zi, vxi, vyi, vzi,
        test_particle_x, test_particle_y, test_particle_z,
        test_particle_vx, test_particle_vy, test_particle_vz,
        test_particle_ax, test_particle_ay, test_particle_az,
        0, # step::Int64
        limit_by_steps # limit_by_steps::Bool
    )
    if is(typeof(sim), Simulation{Cosmological})
        set_zeldovich!(sim, true)
    end
    sim
end

function reset!(s::Simulation)
    # set times and steps
    s.step = 0
    s.t = s.ti
    # set initial positions
    for i in 1:s.w.n
        s.w.particles[i] = withxyz(s.w.particles[i], s.xi[i], s.yi[i], s.zi[i])
    end
    # set initial velocities
    for i in 1:s.w.n
        s.w.vx[i] = s.vxi[i]
        s.w.vy[i] = s.vyi[i]
        s.w.vz[i] = s.vzi[i]
    end
    for i in 1:length(s.test_particle_x)
        s.test_particle_vx[i] = 0.0
        s.test_particle_vy[i] = 0.0
        s.test_particle_vz[i] = 0.0
    end
    if is(typeof(s), Simulation{Cosmological})
        set_zeldovich!(s, true)
    end
end
