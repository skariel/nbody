function save(sim::Simulation, fn::String)
    open(fn, "w") do f
        save(sim, f)
    end
end

function save(sim::Simulation, f::IOStream)
    save(sim.w, f)
    write(f, sim.t)
    write(f, sim.ti)
    write(f, sim.dt)
    write(f, sim.tf)
    write(f, sim.stepf)

    for i in 1:sim.w.n
        write(f, sim.xi[i])
        write(f, sim.yi[i])
        write(f, sim.zi[i])
        write(f, sim.vxi[i])
        write(f, sim.vyi[i])
        write(f, sim.vzi[i])
    end
    tn = length(sim.test_particle_x)
    write(f, tn)

    for i in 1:tn
        write(f, sim.test_particle_x[i])
        write(f, sim.test_particle_y[i])
        write(f, sim.test_particle_z[i])
        write(f, sim.test_particle_vx[i])
        write(f, sim.test_particle_vy[i])
        write(f, sim.test_particle_vz[i])
        write(f, sim.test_particle_ax[i])
        write(f, sim.test_particle_ay[i])
        write(f, sim.test_particle_az[i])
    end

    write(f, sim.step)
    write(f, sim.limit_by_steps)
    write(f, is(typeof(sim.w.space), Cosmological))

    nothing
end

function load_sim(fn::String)
    open(fn, "r") do f
        return load_sim(f)
    end
end

function load_sim(f::IOStream)

    w = load_w(f)
    t = read(f, Float64)
    ti= read(f, Float64)
    dt = read(f, Float64)
    tf = read(f, Float64)
    stepf = read(f, Int64)

    xi = SharedArray(Float64, w.n)
    yi = SharedArray(Float64, w.n)
    zi = SharedArray(Float64, w.n)
    vxi = SharedArray(Float64, w.n)
    vyi = SharedArray(Float64, w.n)
    vzi = SharedArray(Float64, w.n)
    for i in 1:w.n
        xi[i] = read(f, Float64)
        yi[i] = read(f, Float64)
        zi[i] = read(f, Float64)
        vxi[i] = read(f, Float64)
        vyi[i] = read(f, Float64)
        vzi[i] = read(f, Float64)
    end
    tn = read(f, Int64)
    @show tn

    test_particle_x = SharedArray(Float64, tn)
    test_particle_y = SharedArray(Float64, tn)
    test_particle_z = SharedArray(Float64, tn)
    test_particle_vx = SharedArray(Float64, tn)
    test_particle_vy = SharedArray(Float64, tn)
    test_particle_vz = SharedArray(Float64, tn)
    test_particle_ax = SharedArray(Float64, tn)
    test_particle_ay = SharedArray(Float64, tn)
    test_particle_az = SharedArray(Float64, tn)
    for i in 1:tn
        test_particle_x[i] = read(f, Float64)
        test_particle_y[i] = read(f, Float64)
        test_particle_z[i] = read(f, Float64)
        test_particle_vx[i] = read(f, Float64)
        test_particle_vy[i] = read(f, Float64)
        test_particle_vz[i] = read(f, Float64)
        test_particle_ax[i] = read(f, Float64)
        test_particle_ay[i] = read(f, Float64)
        test_particle_az[i] = read(f, Float64)
    end
    step = read(f, Int64)
    limit_by_steps = read(f, Bool)
    is_cosmological = read(f, Bool)

    T = is_cosmological? Cosmological : Newtonian

    sim = Simulation{T}(
        createtree(w), w,
        t, ti, dt, tf, stepf,
        xi, yi, zi,
        vxi, vyi, vzi,
        test_particle_x, test_particle_y, test_particle_z,
        test_particle_vx, test_particle_vy, test_particle_vz,
        test_particle_ax, test_particle_ay, test_particle_az,
        step, limit_by_steps
    )
    return sim
end
