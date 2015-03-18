function World{T<:SpaceType}(particles::Array{Particle, 1}, smth::Number, opening_alpha::Number, dtfrac::Number, space::Type{T})
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
        0.0
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

world{T<:SpaceType}(particles::Array{Particle, 1}; smth=0.01, opening_alpha=0.7, dtfrac=0.035, space::Type{T}=Newtonian) =
    World(particles, smth, opening_alpha, dtfrac, space)

worldnormal{T<:SpaceType}(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0, space::Type{T}=Newtonian) =
    World(particlesnormal(n, scale), smth, opening_alpha, dtfrac, space)

worldspherical{T<:SpaceType}(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0, xy=false, space::Type{T}=Newtonian) =
    World(particlesspherical(n, scale, xy), smth, opening_alpha, dtfrac, space)
