using Base.Test

function calc_accel_brute_force!(w::World, ixs=1:w.n)
    ax,ay,az = calc_accel_brute_force(w, 1:w.n)
    for i in ixs
        w.ax[i] = ax[i]
        w.ay[i] = ay[i]
        w.az[i] = az[i]
    end
end

function calc_accel_brute_force(w::World, ixs=1:w.n)
    ax = zeros(w.n)
    ay = zeros(w.n)
    az = zeros(w.n)
    for i in ixs
        p_i = w.particles[i]
        for j in 1:w.n
            i==j && continue
            pj = w.particles[j]
            const dx = pj._x - p_i._x
            const dy = pj._y - p_i._y
            const dz = pj._z - p_i._z
            const dx2 = dx*dx
            const dy2 = dy*dy
            const dz2 = dz*dz
            const r2 = dx2+dy2+dz2
            const r22 = r2+w.smth2
            const r21 = sqrt(r22)
            const denom = r22*r21/p_i._m
            ax[i] += dx/denom
            ay[i] += dy/denom
            az[i] += dz/denom
        end
    end
    ax[ixs],ay[ixs],az[ixs]
end

function calc_k_energy(w::World, ixs=1:w.n)
    k = 0.0
    @simd for i in 1:w.n
        vx = w.vx[i]
        vy = w.vy[i]
        k += vx*vx+vy*vy
    end
    0.5*w.particles[1]._m*k
end

function calc_p_energy(w::World)
    p = 0.0
    for i in 1:w.n
        p_i = w.particles[i]
        for j in 1:w.n
            i==j && continue
            pj = w.particles[j]
            const dx = pj._x - p_i._x
            const dy = pj._y - p_i._y
            const dz = pj._z - p_i._z
            const dx2 = dx*dx
            const dy2 = dy*dy
            const dz2 = dz*dz
            const r2 = dx2+dy2+dz2
            const denom = sqrt(r2+w.smth2)
            p -= 1.0/denom
        end
    end
    0.5*p*w.particles[1]._m*w.particles[1]._m
end

calc_tot_energy(w::World) = calc_k_energy(w) + calc_p_energy(w)

type TestAcc
    ax_tree::Array{Float64,1}
    a_tree::Array{Float64,1}
    ax_bf::Array{Float64,1}
    a_bf::Array{Float64,1}
    ferr::Array{Float64,1}
    fbelow::Array{Float64,1}
end

function test_acc(n, nout; smth=0.01, opening_alpha=0.7, tp=:spherical, use_brute_force=false)
    ixs = randperm(n)[1:nout]
    w::World
    if tp==:spherical
        w = worldspherical(n, smth=smth, opening_alpha=opening_alpha)
    else
        w = worldnormal(n, smth=smth, opening_alpha=opening_alpha)
    end
    buildtree!(w, createtree(w))

    if !use_brute_force
        for i in ixs
            calculate_accel_on_particle!(w, i)
        end
    else
        calc_accel_brute_force!(w, ixs)
    end

    ax_tree = w.ax[ixs]
    ay_tree = w.ay[ixs]
    az_tree = w.az[ixs]
    a_tree = sqrt(ax_tree.^2+ay_tree.^2+ay_tree.^2)
    println("calculating BF...")
    ax_bf, ay_bf, az_bf = calc_accel_brute_force(w, ixs)
    println("Done!\n")
    a_bf = sqrt(ax_bf.^2+ay_bf.^2+ay_bf.^2)

    dax = ax_tree-ax_bf
    day = ay_tree-ay_bf
    daz = az_tree-az_bf
    ferr = sqrt(dax.^2+day.^2+daz.^2)./a_bf.*100.0;
    sort!(ferr)

    fbelow = [1:nout]./nout.*100.0

    TestAcc(
        ax_tree,
        a_tree,
        ax_bf,
        a_bf,
        ferr,
        fbelow
        )
end

function test_test_particles()
    sim = Simulation(worldnormal(500, smth=0.01, opening_alpha=0.7),
                limit_by_steps=true, stepc=10, n_test_particle=500)

    for i in 1:500
        sim.test_particle_x[i] = sim.w.particles[i]._x
        sim.test_particle_y[i] = sim.w.particles[i]._y
        sim.test_particle_z[i] = sim.w.particles[i]._z
    end
    ox = deepcopy(sim.test_particle_x)
    oy = deepcopy(sim.test_particle_y)
    oz = deepcopy(sim.test_particle_z)

    calc_accel!(sim, Val{true})
    for i in 1:500
        @test sim.test_particle_ax[i] == sim.w.ax[i]
        @test sim.test_particle_ay[i] == sim.w.ay[i]
        @test sim.test_particle_az[i] == sim.w.az[i]
        @test sim.test_particle_ax[i] != 0.0
        @test sim.test_particle_ay[i] != 0.0
        @test sim.test_particle_az[i] != 0.0
    end

    kick!(sim, Val{true}; dt=1.0)
    for i in 1:500
        @test sim.test_particle_vx[i] == sim.w.vx[i]
        @test sim.test_particle_vy[i] == sim.w.vy[i]
        @test sim.test_particle_vz[i] == sim.w.vz[i]
        @test sim.test_particle_vx[i] != 0.0
        @test sim.test_particle_vy[i] != 0.0
        @test sim.test_particle_vz[i] != 0.0
    end

    drift!(sim, Val{true}; dt=1.0)
    for i in 1:500
        @test sim.test_particle_x[i] == sim.w.particles[i]._x
        @test sim.test_particle_y[i] == sim.w.particles[i]._y
        @test sim.test_particle_z[i] == sim.w.particles[i]._z
        @test sim.test_particle_x[i] != 0.0
        @test sim.test_particle_y[i] != 0.0
        @test sim.test_particle_z[i] != 0.0
    end

    for i in 1:500
        sim.test_particle_x[i] = sim.xi[i]
        sim.test_particle_y[i] = sim.yi[i]
        sim.test_particle_z[i] = sim.zi[i]
    end
    exec!(sim, silent=true, test_particles=Val{true})

    for i in 1:500
        @test sim.test_particle_x[i] == sim.w.particles[i]._x
        @test sim.test_particle_y[i] == sim.w.particles[i]._y
        @test sim.test_particle_z[i] == sim.w.particles[i]._z
        @test ox[i] != sim.w.particles[i]._x
        @test oy[i] != sim.w.particles[i]._y
        @test oz[i] != sim.w.particles[i]._z
    end


end


function time_test_particles(n=10000)
    sim = Simulation(worldnormal(n, smth=0.01, opening_alpha=0.7),
                limit_by_steps=true, stepc=100, n_test_particle=3*n)
    for i in 1:n
        sim.test_particle_x[i] = sim.xi[i]
        sim.test_particle_y[i] = sim.yi[i]
        sim.test_particle_z[i] = sim.zi[i]

        sim.test_particle_x[i+n] = sim.xi[i]
        sim.test_particle_y[i+n] = sim.yi[i]
        sim.test_particle_z[i+n] = sim.zi[i]

        sim.test_particle_x[i+2n] = sim.xi[i]
        sim.test_particle_y[i+2n] = sim.yi[i]
        sim.test_particle_z[i+2n] = sim.zi[i]
    end
    t1 = @time exec!(sim; silent=true, test_particles=Val{false})
    for i in 1:n
        sim.test_particle_x[i] = sim.xi[i]
        sim.test_particle_y[i] = sim.yi[i]
        sim.test_particle_z[i] = sim.zi[i]

        sim.test_particle_x[i+n] = sim.xi[i]
        sim.test_particle_y[i+n] = sim.yi[i]
        sim.test_particle_z[i+n] = sim.zi[i]

        sim.test_particle_x[i+2n] = sim.xi[i]
        sim.test_particle_y[i+2n] = sim.yi[i]
        sim.test_particle_z[i+2n] = sim.zi[i]
    end

    t2 = @time exec!(sim, silent=true, test_particles=Val{true})
    @show t2/t1/4
end

function test_grad(opt::Optimization, sim::Simulation; npart=500, nsub=50, stepc=10, scale=1.0)
    ixs = randperm(npart)[1:nsub]

    sim = Simulation(worldspherical(npart, smth=0.1, opening_alpha=0.7, scale=scale),
                limit_by_steps=true, stepc=stepc, n_test_particle=3*npart)

    # getting fast gradients
    opt = Optimization(sim)
    grad!(opt,sim)
    fast_gx = opt.gx[ixs]

    # getting slow gradients
    slow_gx = zeros(nsub)
    exec!(sim; silent=true)
    D = 1.e-8
    g0 = grade(opt, sim)
    for i in 1:nsub
        ix = ixs[i]
        sim.xi[ix] += D
        exec!(sim; silent=true)
        sim.xi[ix] -= D
        g1 = grade(opt, sim)
        slow_gx[i] = (g1-g0)/D
    end

    fast_gx, slow_gx
end

function test_grad(;npart=500, nsub=50, stepc=10, scale=1.0)
    ixs = randperm(npart)[1:nsub]

    sim = Simulation(worldspherical(npart, smth=0.1, opening_alpha=0.7, scale=scale),
                limit_by_steps=true, stepc=stepc, n_test_particle=3*npart)

    # getting fast gradients
    opt = Optimization(sim)
    grad!(opt,sim)
    fast_gx = opt.gx[ixs]

    # getting slow gradients
    slow_gx = zeros(nsub)
    exec!(sim; silent=true)
    D = 1.e-8
    g0 = grade(opt, sim)
    for i in 1:nsub
        ix = ixs[i]
        sim.xi[ix] += D
        exec!(sim; silent=true)
        sim.xi[ix] -= D
        g1 = grade(opt, sim)
        slow_gx[i] = (g1-g0)/D
    end

    fast_gx, slow_gx
end

a=0.3
T=Ta(a)
@test AD(T) == ADa(a)
@test ADD(T) == ADDa(a)
@test H(T) == Ha(a)
@test Z(T) == Za(a);

test_test_particles()
test_grad()

println("test passed!")
