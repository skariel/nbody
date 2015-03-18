function grad!(opt::Optimization)
    # setup test particle positions
    D = 1.e-8
     for i in 1:opt.sim.w.n
        opt.sim.test_particle_x[i] = opt.sim.xi[i] + D
        opt.sim.test_particle_y[i] = opt.sim.yi[i]
        opt.sim.test_particle_z[i] = opt.sim.zi[i]
        opt.sim.test_particle_x[i+opt.sim.w.n] = opt.sim.xi[i]
        opt.sim.test_particle_y[i+opt.sim.w.n] = opt.sim.yi[i] + D
        opt.sim.test_particle_z[i+opt.sim.w.n] = opt.sim.zi[i]
        opt.sim.test_particle_x[i+2opt.sim.w.n] = opt.sim.xi[i]
        opt.sim.test_particle_y[i+2opt.sim.w.n] = opt.sim.yi[i]
        opt.sim.test_particle_z[i+2opt.sim.w.n] = opt.sim.zi[i] + D
    end

    # run the simulation with these test particles
    exec!(opt.sim, true; silent=true)

    # calculate gradient for each particle
     for i in 1:opt.sim.w.n
        dx_o = opt.sim.w.particles[i]._x - opt.x0[i]
        dy_o = opt.sim.w.particles[i]._y - opt.y0[i]
        dz_o = opt.sim.w.particles[i]._z - opt.z0[i]
        d2_o = dx_o*dx_o+dy_o*dy_o+dz_o*dz_o

        dx_n = opt.sim.test_particle_x[i] - opt.x0[i]
        dy_n = opt.sim.test_particle_y[i] - opt.y0[i]
        dz_n = opt.sim.test_particle_z[i] - opt.z0[i]
        opt.gx[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        dx_n = opt.sim.test_particle_x[i+opt.sim.w.n] - opt.x0[i]
        dy_n = opt.sim.test_particle_y[i+opt.sim.w.n] - opt.y0[i]
        dz_n = opt.sim.test_particle_z[i+opt.sim.w.n] - opt.z0[i]
        opt.gy[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        dx_n = opt.sim.test_particle_x[i+2opt.sim.w.n] - opt.x0[i]
        dy_n = opt.sim.test_particle_y[i+2opt.sim.w.n] - opt.y0[i]
        dz_n = opt.sim.test_particle_z[i+2opt.sim.w.n] - opt.z0[i]
        opt.gz[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D
    end
    nothing
end

function grade(opt::Optimization)
    g = 0.0
     for i in 1:opt.sim.w.n
        dx = opt.sim.w.particles[i]._x - opt.x0[i]
        dy = opt.sim.w.particles[i]._y - opt.y0[i]
        dz = opt.sim.w.particles[i]._z - opt.z0[i]
        g += dx*dx+dy*dy+dz*dz
    end
    g
end
