function World(particles::Array{Particle, 1}, smth, opening_excluded_frac, opening_alpha, dtfrac)
    n = length(particles)
    World(
        OctTree(Particle; n=trunc(Integer,4.1*n)),
        particles,
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(n),
        zeros(n),
        n,
        opening_alpha^2,
        opening_excluded_frac^2,
        smth*smth,
        dtfrac
    )
end

particlesnormal(n::Int64) = hilbertsort!([Particle(randn(), randn(), randn(), 1./n) for i in 1:n])

function particlesspherical(n::Int64)
    particles = Particle[]
    while length(particles) < n
        x = 2.0*rand()-1.0
        y = 2.0*rand()-1.0
        z = 2.0*rand()-1.0
        if x*x+y*y+z*z < 1
            push!(particles, Particle(x,y,z,1.0/n))
        end
    end
    hilbertsort!(particles)
end

world(particles::Array{Particle, 1}; smth=0.01, opening_excluded_frac=0.6, opening_alpha=0.7, dtfrac=0.035) =
    World(particles, smth, opening_excluded_frac, opening_alpha, dtfrac)

worldnormal(n::Int64; smth=0.01, opening_excluded_frac=0.6, opening_alpha=0.7, dtfrac=0.035) =
    World(particlesnormal(n), smth, opening_excluded_frac, opening_alpha, dtfrac)

worldspherical(n::Int64; smth=0.01, opening_excluded_frac=0.6, opening_alpha=0.7, dtfrac=0.035) =
    World(particlesspherical(n), smth, opening_excluded_frac, opening_alpha, dtfrac)
