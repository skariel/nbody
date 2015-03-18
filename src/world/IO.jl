
function save(w::World, fn::String)
    open(fn, "w") do f
        save(w, f)
    end
end

function save(w::World, f::IOStream)
    write(f, w.n)
    for p in w.particles
        write(f, p._x)
        write(f, p._y)
        write(f, p._z)
        write(f, p._m)
    end
    for i in 1:w.n
        write(f, w.vx[i])
        write(f, w.vy[i])
        write(f, w.vz[i])
        write(f, w.ax[i])
        write(f, w.ay[i])
        write(f, w.az[i])
        write(f, w.dt[i])
    end
    write(f, w.opening_alpha2)
    write(f, w.smth2)
    write(f, w.mindt)
    write(f, w.dtfrac)
    write(f, w.Ω0)
    write(f, w.ΩΛ)
    write(f, is(typeof(w.space), Cosmological))
    nothing
end

function load_w(fn::String)
    open(fn, "r") do f
        return load_w(f)
    end
end

function load_w(f::IOStream)
    n = read(f, Int64)
    particles = SharedArray(Particle, n)
    for i in 1:n
        particles[i] = Particle(
            read(f, Float64), read(f, Float64),
            read(f, Float64), read(f, Float64)
        )
    end
    vx = SharedArray(Float64, n)
    vy = SharedArray(Float64, n)
    vz = SharedArray(Float64, n)
    ax = SharedArray(Float64, n)
    ay = SharedArray(Float64, n)
    az = SharedArray(Float64, n)
    dt = SharedArray(Float64, n)
    for i in 1:n
        vx[i] = read(f, Float64)
        vy[i] = read(f, Float64)
        vz[i] = read(f, Float64)
        ax[i] = read(f, Float64)
        ay[i] = read(f, Float64)
        az[i] = read(f, Float64)
        dt[i] = read(f, Float64)
    end
    opening_alpha2 = read(f, Float64)
    smth2 = read(f, Float64)
    mindt = read(f, Float64)
    dtfrac = read(f, Float64)

    Ω0 = read(f, Float64)
    ΩΛ = read(f, Float64)


    is_cosmological = read(f, Bool)

    tree = CompiledOctTree(round(Int,1.8*n), Particle)

    T = is_cosmological? Cosmological : Newtonian
    return World{T}(
        tree,
        particles,
        vx, vy, vz,
        ax, ay, az,
        dt, n,
        opening_alpha2,
        smth2, dtfrac,
        T(), mindt, Ω0, ΩΛ
    )
end
