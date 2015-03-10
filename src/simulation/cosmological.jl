function kick!(sim::Simulation{Cosmological}, simulate_test_particles::Bool; dt=0.0)
    # real particles
    a1 = sim.t
    a2 = sim.t+dt
    fk = FK(a1,a2)
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
    a1 = sim.t
    a2 = sim.t+dt
    fd = FD(a1,a2)
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
