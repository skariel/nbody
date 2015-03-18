function calc_dt(sim::Simulation{Cosmological}, simulate_test_particles::Bool)
    mindt = 1.e30 # infinity, ha!
    hat = Ha(sim.t, sim.w)*sim.t
    # real particles
    for i in 1:sim.w.n
        const a2_nonp = sim.w.ax[i]*sim.w.ax[i] + sim.w.ay[i]*sim.w.ay[i] + sim.w.az[i]*sim.w.az[i]
        const a2 = a2_nonp/sim.t/sim.t/sim.t/sim.t # physical accel
        const dyn_dt2 = 2.0*sim.w.dtfrac*sqrt(sim.w.smth2)*sim.t/sqrt(a2)

        const ddt = sqrt(dyn_dt2)*hat
        sim.w.dt[i] = ddt

        if ddt < mindt
            mindt = ddt
        end
    end
    sim.w.mindt=mindt
    # TODO: parametrize this aspect ie maximal timestep in both cosmological and newtonian simulations
    if is(typeof(sim.w), World{Cosmological}) && mindt > 0.05
        mindt = 0.05
    end
    mindt
end



function kick!(sim::Simulation{Cosmological}, simulate_test_particles::Bool; dt=0.0)
    # real particles
    fk = FK(sim.t, sim.t+dt, sim.w)
    for i in 1:sim.w.n
        sim.w.vx[i] += sim.w.ax[i]*fk
        sim.w.vy[i] += sim.w.ay[i]*fk
        sim.w.vz[i] += sim.w.az[i]*fk
    end
    if simulate_test_particles
        # test particles
        for i in 1:length(sim.test_particle_x)
            sim.test_particle_vx[i] += sim.test_particle_ax[i]*fk
            sim.test_particle_vy[i] += sim.test_particle_ay[i]*fk
            sim.test_particle_vz[i] += sim.test_particle_az[i]*fk
        end
    end
    nothing
end

function drift!(sim::Simulation{Cosmological}, simulate_test_particles::Bool; dt=0.0)
    updatespace!(sim.t, sim.w)
    # real particles
    fd = FD(sim.t, sim.t+dt, sim.w)
    for i in 1:sim.w.n
        const dx = sim.w.vx[i]*fd
        const dy = sim.w.vy[i]*fd
        const dz = sim.w.vz[i]*fd
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    if simulate_test_particles
        # test particles
        for i in 1:length(sim.test_particle_x)
            sim.test_particle_x[i] += sim.test_particle_vx[i]*fd
            sim.test_particle_y[i] += sim.test_particle_vy[i]*fd
            sim.test_particle_z[i] += sim.test_particle_vz[i]*fd
        end
    end
    nothing
end

updatespace!(t::Float64, w::World{Cosmological}) = nothing

function calc_accel!(sim::Simulation{Cosmological}, simulate_test_particles::Bool)
    updatespace!(sim.t, sim.w)
    buildtree!(sim.w, sim.tree)
    if simulate_test_particles
        calc_accel!(sim.w,
            sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
            sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az)
    else
        calc_accel!(sim.w)
    end
    nothing
end

function exec!(sim::Simulation{Cosmological}, simulate_test_particles::Bool; use_brute_force=false, silent=false)
    reset!(sim)
    tic()
    calc_accel!(sim, simulate_test_particles)
    break_time = false
    while true
        sim.step += 1
        if sim.step % 500 == 0
            @show sim.dt
            @show sim.step
            @show sim.t
        end
        # physical dt!
        sim.dt = calc_dt(sim, simulate_test_particles)
        if !sim.limit_by_steps
            if sim.t+sim.dt > sim.tf
                sim.dt = sim.tf - sim.t
                break_time = true
            end
        end
        !silent && print("s=",sim.step," t=",sim.t," dt=",sim.dt)
        tic()

        kick!(sim, simulate_test_particles, dt=sim.dt/2)
        drift!(sim, simulate_test_particles, dt=sim.dt)
        if !use_brute_force
            calc_accel!(sim, simulate_test_particles)
        else
            calc_accel_brute_force!(sim)
        end
        sim.t += sim.dt/2
        kick!(sim, simulate_test_particles, dt=sim.dt/2)
        sim.t += sim.dt/2

        elapsed = toq()
        !silent && print(" /",elapsed,"s\n")

        sim.limit_by_steps && sim.step >= sim.stepf && break
        break_time && break
    end
    !silent && println("\n--- Done!\n")
    silent? toq() : toc()
    nothing
end

function set_zeldovich!(sim::Simulation, simulate_test_particles::Bool)
    calc_accel!(sim, simulate_test_particles)
    ah = sim.t*Ha(sim.t, sim.w)
    fac1 = 2.0/3.0/ah
    for i in 1:sim.w.n
        sim.w.vxi[i] = sim.w.ax[i]*fac1
        sim.w.vyi[i] = sim.w.ay[i]*fac1
        sim.w.vzi[i] = sim.w.az[i]*fac1
    end
    if simulate_test_particles
        sim.test_particle_vx[i] = sim.test_particle_ax[i]*fac1
        sim.test_particle_vy[i] = sim.test_particle_ay[i]*fac1
        sim.test_particle_vz[i] = sim.test_particle_az[i]*fac1
    end
end
