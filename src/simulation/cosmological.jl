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
        s.w.vx[i] = s.w.ax[i] / s.w.space.F3
        s.w.vy[i] = s.w.ay[i] / s.w.space.F3
        s.w.vz[i] = s.w.az[i] / s.w.space.F3
    end
    for i in 1:length(s.test_particle_x)
        s.test_particle_vx[i] = s.test_particle_ax[i] / s.w.space.F3
        s.test_particle_vy[i] = s.test_particle_ay[i] / s.w.space.F3
        s.test_particle_vz[i] = s.test_particle_az[i] / s.w.space.F3
    end
end

function drift!(sim::Simulation{Cosmological}; dt=0.0)
    updatespace!(sim.t, sim.w)
    # real particles
    @inline for i in 1:sim.w.n
        const dx = sim.w.vx[i]*dt*sim.w.space.sa
        const dy = sim.w.vy[i]*dt*sim.w.space.sa
        const dz = sim.w.vz[i]*dt*sim.w.space.sa
        sim.w.particles[i] = addxyz(sim.w.particles[i], dx, dy, dz)
    end
    # test particles
    @inline for i in 1:length(sim.test_particle_x)
        sim.test_particle_x[i] += sim.test_particle_vx[i]*dt*sim.w.space.sa
        sim.test_particle_y[i] += sim.test_particle_vy[i]*dt*sim.w.space.sa
        sim.test_particle_z[i] += sim.test_particle_vz[i]*dt*sim.w.space.sa
    end
    nothing
end

@inline updatespace!(t::Float64, w::World{Cosmological}) = (w.space=Cosmological(t); nothing)

function calc_accel!(sim::Simulation{Cosmological})
    updatespace!(sim.t, sim.w)
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w,
        sim.test_particle_x, sim.test_particle_y, sim.test_particle_z,
        sim.test_particle_ax, sim.test_particle_ay, sim.test_particle_az,
        sim.test_particle_vx, sim.test_particle_vy, sim.test_particle_vz)
    nothing
end
