
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




function drift!(sim::Simulation{Newtonian}; dt=0.0)
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

function calc_accel!(sim::Simulation{Newtonian})
    updatespace!(sim.t, sim.w)
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w,
        sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
        sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az)
    nothing
end
