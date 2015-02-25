function exec(p::Particle, hist::History)
    vx = 0.0
    vy = 0.0
    vz = 0.0
    ax = 0.0
    ay = 0.0
    az = 0.0
    tree = hist.tree[1]
    ax, ay, az = calc_accel(p,tree,hist.w)
    for i in 1:(length(hist.dt)-1)
        # kick half dt
        dt = hist.dt[i]
        vx += ax*dt/2
        vy += ay*dt/2
        vz += az*dt/2
        # drift dt
        p = addxyz(p, vx*dt, vy*dt, vz*dt)
        # kick another half dt
        tree = hist.tree[i+1]
        ax, ay, az = calc_accel(p,tree,hist.w)
        vx += ax*dt/2
        vy += ay*dt/2
        vz += az*dt/2
    end
    p
end

@inline function grade(pf::Particle, p0::Particle)
    const dx = pf._x - p0._x
    const dy = pf._y - p0._y
    const dz = pf._z - p0._z
    dx*dx+dy*dy+dz*dz
end

function grade(sim::Simulation, opt::Optimization)
    g = 0.0
    @inbounds for i in 1:sim.w.n
        g += grade(sim.w.pf[i], Particle(opt.x0[i], opt.y0[i], opy.z0[i]))
    end
    g
end

@inline function grad(pi::Particle, p0::Particle, pf::Particle, hist::History)
    const D = 1.e-6

    const new_p_x = exec(addxyz(pi,  D,0.0,0.0), opt.hist)
    const new_p_y = exec(addxyz(pi,0.0,  D,0.0), opt.hist)
    const new_p_z = exec(addxyz(pi,0.0,0.0,  D), opt.hist)

    const org_grade = grade(pf, p0)
    const x_grade = grade(new_p_x, p0)
    const y_grade = grade(new_p_y, p0)
    const z_grade = grade(new_p_z, p0)

    (x_grade-org_grade)/D, (y_grade-org_grade)/D, (z_grade-org_grade)/D
end

function grad!(rng::UnitRange{Int64}, opt::Optimization, sim::Simulation)
    for i in rng
        @inbounds pi = Particle(sim.xi[i]. sim.yi[i], sim.zi[i])
        @inbounds p0 = Particle(opt.x0[i]. opt.y0[i], opt.z0[i])
        @inbounds pf = sim.w.pf[i]
        @inbounds opt.gx[i], opt.gy[i], opt.gz[i] = grad(pi, p0, pf, hist)
    end
end

function grad!(opt::Optimization, sim::Simulation)
    chunks = get_chunks(sim.w.n)
    @sync for i in 1:length(workers())
        @async remotecall_wait(workers()[i], grad!, chunks[i], opt, sim)
    end    
end
