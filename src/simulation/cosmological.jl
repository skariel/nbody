function calc_dt(sim::Simulation{Cosmological}, simulate_test_particles::Bool)
    mindt2 = 1.e30 # infinity, ha!
    hat = Ha(sim.t)*sim.t
    # real particles
    @inline for i in 1:sim.w.n
        const a2_nonp = sim.w.ax[i]*sim.w.ax[i] + sim.w.ay[i]*sim.w.ay[i] + sim.w.az[i]*sim.w.az[i]
        const a2 = a2_nonp/sim.t/sim.t/sim.t/sim.t # physical accel
        const dyn_dt2 = 2.0*sim.w.dtfrac*sqrt(sim.w.smth2)*sim.t/sqrt(a2)
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    if simulate_test_particles
        # test particles
        @inline for i in 1:length(sim.test_particle_x)
            const a2_nonp = sim.test_particle_ax[i]*sim.test_particle_ax[i] + sim.test_particle_ay[i]*sim.test_particle_ay[i] + sim.test_particle_az[i]*sim.test_particle_az[i]
            const dyn_dt2 = sqrt(sim.w.smth2/a2)*sim.w.dtfrac*sim.w.dtfrac
            const a2 = a2_nonp/sim.t/sim.t/sim.t/sim.t # physical accel
            const dyn_dt2 = 2.0*sim.w.dtfrac*sqrt(sim.w.smth2)*sim.t/sqrt(a2)
            if dyn_dt2 < mindt2
                mindt2 = dyn_dt2
            end
        end
    end
    mindt = sqrt(mindt2)*hat
    # TODO: parametrize this aspect ie maximal timestep in both cosmological and newtonian simulations
    if is(typeof(sim.w), World{Cosmological}) && mindt > 0.05
        mindt = 0.05
    end
    mindt
end



function kick!(sim::Simulation{Cosmological}, simulate_test_particles::Bool; dt=0.0)
    # real particles
    fk = FK(sim.t, sim.t+dt)
    @inline for i in 1:sim.w.n
        sim.w.vx[i] += sim.w.ax[i]*fk
        sim.w.vy[i] += sim.w.ay[i]*fk
        sim.w.vz[i] += sim.w.az[i]*fk
    end
    if simulate_test_particles
        # test particles
        @inline for i in 1:length(sim.test_particle_x)
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
    fd = FD(sim.t, sim.t+dt)
    @inline for i in 1:sim.w.n
        const dx = sim.w.vx[i]*fd
        const dy = sim.w.vy[i]*fd
        const dz = sim.w.vz[i]*fd
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    if simulate_test_particles
        # test particles
        @inline for i in 1:length(sim.test_particle_x)
            sim.test_particle_x[i] += sim.test_particle_vx[i]*fd
            sim.test_particle_y[i] += sim.test_particle_vy[i]*fd
            sim.test_particle_z[i] += sim.test_particle_vz[i]*fd
        end
    end
    nothing
end

@inline updatespace!(t::Float64, w::World{Cosmological}) = (w.space=Cosmological(t); nothing)

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
