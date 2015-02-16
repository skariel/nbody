function worldnormal(n::Int64; smth=0.000001, opening_excluded_frac=0.0, opening_alpha=0.7, dtfrac=0.035)
    particles = [Particle(randn(), randn(), randn(), 1./n) for i in 1:n]
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

function worldspherical(n::Int64; smth=0.0, opening_excluded_frac=0.6, opening_alpha=0.7, dtfrac=0.035)
    particles = Particle[]
    while length(particles) < n
        tx = rand()*2.0-1.0
        ty = rand()*2.0-1.0
        tz = rand()*2.0-1.0
        tx*tx+ty*ty+tz*tz < 1.0 && push!(particles, Particle(tx,ty,tz, 1./n))
    end
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
