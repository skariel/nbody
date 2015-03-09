type Simulation{T<:SpaceType}
    tree::OctTree{Particle}
    w::World{T}
    t::Float64 # time parameter: real time for newtonian space, scale factor `a` for cosmological space
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
        @inbounds xi[i] = w.particles[i]._x
        @inbounds yi[i] = w.particles[i]._y
        @inbounds zi[i] = w.particles[i]._z
        @inbounds vxi[i] = 0.0
        @inbounds vyi[i] = 0.0
        @inbounds vzi[i] = 0.0
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

    Simulation{T}(
        createtree(w),
        w,
        ti,  # t
        0.0, # dt
        ti,  # ti
        tf,  # tf
        stepc,   # stepf::Int64
        xi, yi, zi, vxi, vyi, vzi,
        test_particle_x, test_particle_y, test_particle_z,
        test_particle_vx, test_particle_vy, test_particle_vz,
        test_particle_ax, test_particle_ay, test_particle_az,
        0, # step::Int64
        limit_by_steps # limit_by_steps::Bool
    )
end

function reset_vel!(s::Simulation{Newtonian})
    for i in 1:s.w.n
        s.w.vx[i] = 0.0
        s.w.vy[i] = 0.0
        s.w.vz[i] = 0.0
        s.vxi[i] = 0.0
        s.vyi[i] = 0.0
        s.vzi[i] = 0.0
    end
    for i in 1:length(s.test_particle_x)
        s.test_particle_vx[i] = 0.0
        s.test_particle_vy[i] = 0.0
        s.test_particle_vz[i] = 0.0
    end
end

function reset_vel!(s::Simulation{Cosmological})
    # use Zeldovich approx:

    # 1) zero velocities
    for i in 1:s.w.n
        s.w.vx[i] = 0.0
        s.w.vy[i] = 0.0
        s.w.vz[i] = 0.0
        s.vxi[i] = 0.0
        s.vyi[i] = 0.0
        s.vzi[i] = 0.0
    end
    for i in 1:length(s.test_particle_x)
        s.test_particle_vx[i] = 0.0
        s.test_particle_vy[i] = 0.0
        s.test_particle_vz[i] = 0.0
    end

    # 2) calculate accel
    calc_accel!(s)

    # 3) fix Zeldovich velocities
    for i in 1:s.w.n
        s.w.vx[i] = s.w.ax[i] / s.w.space.FU
        s.w.vy[i] = s.w.ay[i] / s.w.space.FU
        s.w.vz[i] = s.w.az[i] / s.w.space.FU
    end
    for i in 1:length(s.test_particle_x)
        s.test_particle_vx[i] = s.test_particle_ax[i] / s.w.space.FU
        s.test_particle_vy[i] = s.test_particle_ay[i] / s.w.space.FU
        s.test_particle_vz[i] = s.test_particle_az[i] / s.w.space.FU
    end
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
    reset_vel!(s)
end

function calc_dt(sim::Simulation)
    mindt2 = 1.e30 # infinity, ha!
    # real particles
    @inline for i in 1:sim.w.n
        const a2 = sim.w.ax[i]*sim.w.ax[i] + sim.w.ay[i]*sim.w.ay[i] + sim.w.az[i]*sim.w.az[i]
        const dyn_dt2 = sqrt(sim.w.smth2/a2)*sim.w.dtfrac*sim.w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    # test particles
    @inline for i in 1:length(sim.test_particle_x)
        const a2 = sim.test_particle_ax[i]*sim.test_particle_ax[i] + sim.test_particle_ay[i]*sim.test_particle_ay[i] + sim.test_particle_az[i]*sim.test_particle_az[i]
        const dyn_dt2 = sqrt(sim.w.smth2/a2)*sim.w.dtfrac*sim.w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    mindt = sqrt(mindt2)
    # TODO: parametrize this aspect ie maximal timestep in both cosmological and newtonian simulations
    if is(typeof(sim.w), World{Cosmological}) && mindt > 0.05
        mindt = 0.05
    end
    mindt
end

function kick!(sim::Simulation; dt=0.0)
    # real particles
    @inline for i in 1:sim.w.n
        sim.w.vx[i] += sim.w.ax[i]*dt
        sim.w.vy[i] += sim.w.ay[i]*dt
        sim.w.vz[i] += sim.w.az[i]*dt
    end
    # test particles
    @inline for i in 1:length(sim.test_particle_x)
        sim.test_particle_vx[i] += sim.test_particle_ax[i]*dt
        sim.test_particle_vy[i] += sim.test_particle_ay[i]*dt
        sim.test_particle_vz[i] += sim.test_particle_az[i]*dt
    end
    nothing
end

function drift!(sim::Simulation; dt=0.0)
    # real particles
    @inline for i in 1:sim.w.n
        const dx = sim.w.vx[i]*dt
        const dy = sim.w.vy[i]*dt
        const dz = sim.w.vz[i]*dt
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    # test particles
    @inline for i in 1:length(sim.test_particle_x)
        sim.test_particle_x[i] += sim.test_particle_vx[i]*dt
        sim.test_particle_y[i] += sim.test_particle_vy[i]*dt
        sim.test_particle_z[i] += sim.test_particle_vz[i]*dt
    end
    nothing
end

@inline updatespace!(t::Float64, w::World{Newtonian}) = nothing
@inline updatespace!(t::Float64, w::World{Cosmological}) = (w.space=Cosmological(t); nothing)

function calc_accel!(sim::Simulation{Newtonian})
    updatespace!(sim.t, sim.w)
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w,
        sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
        sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az)
    nothing
end

function calc_accel!(sim::Simulation{Cosmological})
    updatespace!(sim.t, sim.w)
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w,
        sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
        sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az,
        sim.test_particle_vx, sim.test_particle_vy, sim.test_particle_vz)
    nothing
end

function exec!(sim::Simulation; use_brute_force=false, silent=false)
    reset!(sim)
    tic()
    calc_accel!(sim)
    break_time = false
    while true
        sim.step += 1
        sim.dt = calc_dt(sim)
        if !sim.limit_by_steps
            if sim.t+sim.dt > sim.tf
                sim.dt = sim.tf - sim.t
                break_time = true
            end
        end
        !silent && print("s=",sim.step," t=",sim.t," dt=",sim.dt)
        tic()

        kick!(sim, dt=sim.dt/2)
        drift!(sim, dt=sim.dt)
        if !use_brute_force
            calc_accel!(sim)
        else
            calc_accel_brute_force!(sim)
        end
        kick!(sim, dt=sim.dt/2)

        sim.t += sim.dt
        elapsed = toq()
        !silent && print(" /",elapsed,"s\n")

        sim.limit_by_steps && sim.step >= sim.stepf && break
        break_time && break
    end
    !silent && println("\n--- Done!\n")
    silent? toq() : toc()
    nothing
end
