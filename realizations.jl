function World(particles::Array{Particle, 1}, smth, opening_alpha, dtfrac)
    n = length(particles)
    spar = SharedArray(Particle, n)
    vx = SharedArray(Float64, n)
    vy = SharedArray(Float64, n)
    vz = SharedArray(Float64, n)
    ax = SharedArray(Float64, n)
    ay = SharedArray(Float64, n)
    az = SharedArray(Float64, n)

    for i in 1:n
        @inbounds spar[i] = particles[i]
        @inbounds vx[i] = 0.0
        @inbounds vy[i] = 0.0
        @inbounds vz[i] = 0.0
        @inbounds ax[i] = 0.0
        @inbounds ay[i] = 0.0
        @inbounds az[i] = 0.0
    end
    World(
        CompiledOctTree(int(1.8*n), Particle),
        spar,
        vx,vy,vz,ax,ay,az,
        n,
        opening_alpha^2,
        smth*smth,
        dtfrac
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

world(particles::Array{Particle, 1}; smth=0.01, opening_alpha=0.7, dtfrac=0.035) =
    World(particles, smth, opening_alpha, dtfrac)

worldnormal(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0) =
    World(particlesnormal(n, scale), smth, opening_alpha, dtfrac)

worldspherical(n::Int64; smth=0.01, opening_alpha=0.7, dtfrac=0.035, scale=1.0, xy=false) =
    World(particlesspherical(n, scale, xy), smth, opening_alpha, dtfrac)
