function calc_dt(sim::Simulation, ::Type{Val{false}})
    mindt2 = 1.e30 # infinity, ha!
    @inline for i in 1:sim.w.n
        const a2 = sim.w.ax[i]*sim.w.ax[i] + sim.w.ay[i]*sim.w.ay[i] + sim.w.az[i]*sim.w.az[i]
        const dyn_dt2 = sqrt(sim.w.smth2/a2)*sim.w.dtfrac*sim.w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    sqrt(mindt2)
end

function calc_dt(sim::Simulation, ::Type{Val{true}})
    mindt = calc_dt(sim, Val{false})
    mindt2 = mindt*mindt
    @inline for i in 1:length(sim.test_particle_x)
        const a2 = sim.test_particle_ax[i]*sim.test_particle_ax[i] + sim.test_particle_ay[i]*sim.test_particle_ay[i] + sim.test_particle_az[i]*sim.test_particle_az[i]
        const dyn_dt2 = sqrt(sim.w.smth2/a2)*sim.w.dtfrac*sim.w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    sqrt(mindt2)
end

function kick!(sim::Simulation, ::Type{Val{false}}; dt=0.0)
    @inline for i in 1:sim.w.n
        sim.w.vx[i] += sim.w.ax[i]*dt
        sim.w.vy[i] += sim.w.ay[i]*dt
        sim.w.vz[i] += sim.w.az[i]*dt
    end
    nothing
end

function kick!(sim::Simulation, ::Type{Val{true}}; dt=0.0)
    kick!(sim, Val{false}; dt=dt)
    @inline for i in 1:length(sim.test_particle_x)
        sim.test_particle_vx[i] += sim.test_particle_ax[i]*dt
        sim.test_particle_vy[i] += sim.test_particle_ay[i]*dt
        sim.test_particle_vz[i] += sim.test_particle_az[i]*dt
    end
    nothing
end

function drift!(sim::Simulation, ::Type{Val{false}}; dt=0.0)
    @inline for i in 1:sim.w.n
        const dx = sim.w.vx[i]*dt
        const dy = sim.w.vy[i]*dt
        const dz = sim.w.vz[i]*dt
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    nothing
end

function drift!(sim::Simulation, ::Type{Val{true}}; dt=0.0)
    drift!(sim, Val{false}; dt=dt)
    @inline for i in 1:length(sim.test_particle_x)
        sim.test_particle_x[i] += sim.test_particle_vx[i]*dt
        sim.test_particle_y[i] += sim.test_particle_vy[i]*dt
        sim.test_particle_z[i] += sim.test_particle_vz[i]*dt
    end
    nothing
end

function calc_accel!(sim::Simulation, ::Type{Val{false}})
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w)
    nothing
end

function calc_accel!(sim::Simulation, ::Type{Val{true}})
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w,
        sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
        sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az)
    nothing
end

function exec!(sim::Simulation; use_brute_force=false, silent=false, test_particles=Val{false})
    reset!(sim)
    tic()
    calc_accel!(sim, test_particles)
    break_time = false
    while true
        sim.step += 1
        sim.dt = calc_dt(sim, test_particles)
        if !sim.limit_by_steps
            if sim.t+sim.dt > sim.tf
                sim.dt = sim.tf - sim.t
                break_time = true
            end
        end
        !silent && print("s=",sim.step," t=",sim.t," dt=",sim.dt)
        tic()

        kick!(sim, test_particles, dt=sim.dt/2)
        drift!(sim, test_particles, dt=sim.dt)
        if !use_brute_force
            calc_accel!(sim, test_particles)
        else
            calc_accel_brute_force!(sim, test_particles)
        end
        kick!(sim, test_particles, dt=sim.dt/2)

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
