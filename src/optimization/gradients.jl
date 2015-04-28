function grad!(opt::Optimization)
    # setup test particle positions
    D = sqrt(opt.sim.w.smth2)/200000
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
        zx, zy, zz = get_z_coords(
            opt.sim.w.particles[i]._x,
            opt.sim.w.particles[i]._y,
            opt.sim.w.particles[i]._z,
            opt.sim.w.vx[i], opt.sim.w.vy[i], opt.sim.w.vz[i]
        )
        dx_o = zx - opt.x0[i]
        dy_o = zy - opt.y0[i]
        dz_o = zz - opt.z0[i]
        d2_o = dx_o*dx_o+dy_o*dy_o+dz_o*dz_o

        zx, zy, zz = get_z_coords(
            opt.sim.test_particle_x[i],
            opt.sim.test_particle_y[i],
            opt.sim.test_particle_z[i],
            opt.sim.test_particle_vx[i],
            opt.sim.test_particle_vy[i],
            opt.sim.test_particle_vz[i]
        )
        dx_n = zx - opt.x0[i]
        dy_n = zy - opt.y0[i]
        dz_n = zz - opt.z0[i]
        opt.gx[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        zx, zy, zz = get_z_coords(
            opt.sim.test_particle_x[i+opt.sim.w.n],
            opt.sim.test_particle_y[i+opt.sim.w.n],
            opt.sim.test_particle_z[i+opt.sim.w.n],
            opt.sim.test_particle_vx[i+opt.sim.w.n],
            opt.sim.test_particle_vy[i+opt.sim.w.n],
            opt.sim.test_particle_vz[i+opt.sim.w.n]
        )
        dx_n = zx - opt.x0[i]
        dy_n = zy - opt.y0[i]
        dz_n = zz - opt.z0[i]
        opt.gy[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D

        zx, zy, zz = get_z_coords(
            opt.sim.test_particle_x[i+2opt.sim.w.n],
            opt.sim.test_particle_y[i+2opt.sim.w.n],
            opt.sim.test_particle_z[i+2opt.sim.w.n],
            opt.sim.test_particle_vx[i+2opt.sim.w.n],
            opt.sim.test_particle_vy[i+2opt.sim.w.n],
            opt.sim.test_particle_vz[i+2opt.sim.w.n]
        )
        dx_n = zx - opt.x0[i]
        dy_n = zy - opt.y0[i]
        dz_n = zz - opt.z0[i]
        opt.gz[i] = (dx_n*dx_n+dy_n*dy_n+dz_n*dz_n - d2_o)/D
    end
    nothing
end

function grade(opt::Optimization)
    for i in 1:opt.sim.w.n
        zx, zy, zz = get_z_coords(
            opt.sim.w.particles[i]._x,
            opt.sim.w.particles[i]._y,
            opt.sim.w.particles[i]._z,
            opt.sim.w.vx[i], opt.sim.w.vy[i], opt.sim.w.vz[i])
        dx = zx - opt.x0[i]
        dy = zy - opt.y0[i]
        dz = zz - opt.z0[i]
        opt.grad[i] = dx*dx+dy*dy+dz*dz
    end
    mean(opt.grad)
end
