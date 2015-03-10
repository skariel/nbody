function kick!(sim::Simulation{Newtonian}, simulate_test_particles::Bool; dt=0.0)
    # real particles
    @inline for i in 1:sim.w.n
        sim.w.vx[i] += sim.w.ax[i]*dt
        sim.w.vy[i] += sim.w.ay[i]*dt
        sim.w.vz[i] += sim.w.az[i]*dt
    end
    if simulate_test_particles
        # test particles
        @inline for i in 1:length(sim.test_particle_x)
            sim.test_particle_vx[i] += sim.test_particle_ax[i]*dt
            sim.test_particle_vy[i] += sim.test_particle_ay[i]*dt
            sim.test_particle_vz[i] += sim.test_particle_az[i]*dt
        end
    end
    nothing
end


function drift!(sim::Simulation{Newtonian}, simulate_test_particles::Bool; dt=0.0)
    # real particles
    @inline for i in 1:sim.w.n
        const dx = sim.w.vx[i]*dt
        const dy = sim.w.vy[i]*dt
        const dz = sim.w.vz[i]*dt
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    if simulate_test_particles
        # test particles
        @inline for i in 1:length(sim.test_particle_x)
            sim.test_particle_x[i] += sim.test_particle_vx[i]*dt
            sim.test_particle_y[i] += sim.test_particle_vy[i]*dt
            sim.test_particle_z[i] += sim.test_particle_vz[i]*dt
        end
    end
    nothing
end


@inline updatespace!(t::Float64, w::World{Newtonian}) = nothing

function calc_accel!(sim::Simulation{Newtonian}, simulate_test_particles::Bool)
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
