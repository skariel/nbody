function grad!(opt::Optimization, sim::Simulation)
    # setup test particle positions
    D = 1.e-7
    @inbounds for i in 1:sim.w.n
        sim.test_particle_x[i] = sim.xi[i] + D
        sim.test_particle_y[i] = sim.yi[i]
        sim.test_particle_z[i] = sim.zi[i]
        sim.test_particle_x[i+sim.w.n] = sim.xi[i]
        sim.test_particle_y[i+sim.w.n] = sim.yi[i] + D
        sim.test_particle_z[i+sim.w.n] = sim.zi[i]
        sim.test_particle_x[i+2sim.w.n] = sim.xi[i]
        sim.test_particle_y[i+2sim.w.n] = sim.yi[i]
        sim.test_particle_z[i+2sim.w.n] = sim.zi[i] + D
    end

    # run the simulation with these test particles
    exec!(sim; silent=true, test_particles=Val{true})

    # calculate gradient for each particle
    @inbounds for i in 1:sim.w.n
        dx_o = sim.w.particles[i]._x - opt.x0[i]
        dy_o = sim.w.particles[i]._y - opt.y0[i]
        dz_o = sim.w.particles[i]._z - opt.z0[i]
        d2_o = dx_o*dx_o+dy_o*dy_o+dz_o*dz_o

        dx_n = sim.test_particle_x[i] - opt.x0[i]
        dy_n = sim.test_particle_y[i] - opt.y0[i]
        dz_n = sim.test_particle_z[i] - opt.z0[i]
        opt.gx[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        dx_n = sim.test_particle_x[i+sim.w.n] - opt.x0[i]
        dy_n = sim.test_particle_y[i+sim.w.n] - opt.y0[i]
        dz_n = sim.test_particle_z[i+sim.w.n] - opt.z0[i]
        opt.gy[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        dx_n = sim.test_particle_x[i+2sim.w.n] - opt.x0[i]
        dy_n = sim.test_particle_y[i+2sim.w.n] - opt.y0[i]
        dz_n = sim.test_particle_z[i+2sim.w.n] - opt.z0[i]
        opt.gz[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D
    end
    nothing
end

function grade(opt::Optimization, sim::Simulation)
    g = 0.0
    @inbounds for i in 1:sim.w.n
        dx = sim.w.particles[i]._x - opt.x0[i]
        dy = sim.w.particles[i]._y - opt.y0[i]
        dz = sim.w.particles[i]._z - opt.z0[i]
        g += dx*dx+dy*dy+dz*dz
    end
    g
end
