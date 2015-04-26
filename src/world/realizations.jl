function World{T<:SpaceType}(particles::Array{Particle, 1}, smth::Number, opening_alpha::Number, dtfrac::Number, space::Type{T}, Ω0::Float64, ΩΛ::Float64)
    n = length(particles)
    spar = SharedArray(Particle, n)
    vx = SharedArray(Float64, n)
    vy = SharedArray(Float64, n)
    vz = SharedArray(Float64, n)
    ax = SharedArray(Float64, n)
    ay = SharedArray(Float64, n)
    az = SharedArray(Float64, n)
    dt = SharedArray(Float64, n)

    for i in 1:n
         spar[i] = particles[i]
         vx[i] = 0.0
         vy[i] = 0.0
         vz[i] = 0.0
         ax[i] = 0.0
         ay[i] = 0.0
         az[i] = 0.0
         dt[i] = 0.0
    end
    World{space}(
        CompiledOctTree(round(Int,1.8*n), Particle),
        spar,
        vx,vy,vz,ax,ay,az, dt,
        n,
        opening_alpha^2,
        smth*smth,
        dtfrac,
        T(),
        0.0,
        Ω0,
        ΩΛ
    )
end

particlesnormal(n::Int64, scale=1.0) = hilbertsort!([Particle(randn()*scale, randn()*scale, randn()*scale, 1./n) for i in 1:n])

function particlesspherical(n::Int64, scale=1.0, xy=false)
    particles = Particle[]
    while length(particles) < n
        x = 2.0*rand()*scale-scale
        y = 2.0*rand()*scale-scale
        z = xy ? 0.0 : 2.0*rand()*scale-scale
        if x*x+y*y+z*z < scale*scale
            push!(particles, Particle(x,y,z,1.0/n))
        end
    end
    hilbertsort!(particles)
end

world{T<:SpaceType}(particles::Array{Particle, 1}; smth=0.01, opening_alpha=0.7, dtfrac=0.035, space::Type{T}=Newtonian, Ω0=1.0, ΩΛ=0.0) =
    World(particles, smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)

worldnormal{T<:SpaceType}(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0, space::Type{T}=Newtonian, Ω0=1.0, ΩΛ=0.0) =
    World(particlesnormal(n, scale), smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)

worldspherical{T<:SpaceType}(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0, xy=false, space::Type{T}=Newtonian, Ω0=1.0, ΩΛ=0.0) =
    World(particlesspherical(n, scale, xy), smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)

function filter_r{T<:SpaceType}(w::World{T}, r::Float64)
    particle_arr = Particle[]
    vx = Float64[]
    vy = Float64[]
    vz = Float64[]
    rt_2 = r*r
    xs = [p._x for p in w.particles]
    d = 0.5*(maximum(xs) + minimum(xs))
    for i in 1:w.n
        p = w.particles[i]
        x = p._x - d
        y = p._y - d
        z = p._z - d
        r2 = x*x+y*y+z*z
        if r2 < rt_2
            push!(particle_arr, Particle(x, y, z, p._m))
            push!(vx, w.vx[i])
            push!(vy, w.vy[i])
            push!(vz, w.vz[i])
        end
    end
    mvx = mean(vx)
    mvy = mean(vy)
    mvz = mean(vz)
    @show smth = sqrt(w.smth2)
    @show opening_alpha = sqrt(w.opening_alpha2)
    @show dtfrac = w.dtfrac
    @show space = typeof(w.space)
    @show Ω0 = w.Ω0
    @show ΩΛ = w.ΩΛ
    new_world = World(particle_arr, smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)
    for i in 1:new_world.n
        new_world.vx[i] = vx[i]-mvx
        new_world.vy[i] = vy[i]-mvy
        new_world.vz[i] = vz[i]-mvz
    end
    new_world
end

function filter_rate{T<:SpaceType}(w::World{T}, rate::Int64)
    particle_arr = Particle[]
    vx = Float64[]
    vy = Float64[]
    vz = Float64[]
    for i in 1:rate:w.n
        p = w.particles[i]
        push!(particle_arr, p)
        push!(vx, w.vx[i])
        push!(vy, w.vy[i])
        push!(vz, w.vz[i])
    end
    # fixing mass
    for i in 1:length(particle_arr)
        p = particle_arr[i]
        particle_arr[i] = Particle(p._x,p._y,p._z,p._m*w.n/length(particle_arr))
    end
    mvx = mean(vx)
    mvy = mean(vy)
    mvz = mean(vz)
    @show smth = sqrt(w.smth2)
    @show opening_alpha = sqrt(w.opening_alpha2)
    @show dtfrac = w.dtfrac
    @show space = typeof(w.space)
    @show Ω0 = w.Ω0
    @show ΩΛ = w.ΩΛ
    new_world = World(particle_arr, smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)
    for i in 1:new_world.n
        new_world.vx[i] = vx[i]-mvx
        new_world.vy[i] = vy[i]-mvy
        new_world.vz[i] = vz[i]-mvz
    end
    new_world
end
