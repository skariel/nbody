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

function test_backdynamics()
    sim = Simulation(worldnormal(500, smth=0.01, opening_alpha=0.7),
                limit_by_steps=true, stepc=10)

    op_i = sim.w.particles[1]

    hist = exec!(sim, silent=true, accumulate_history=true)

    op_f = sim.w.particles[1]

    odx = op_f._x - op_i._x
    ody = op_f._y - op_i._y
    odz = op_f._z - op_i._z
    org_d = sqrt(odx*odx+ody*ody+odz*odz)

    np_f = exec(op_i, hist)
    ndx = np_f._x - op_i._x
    ndy = np_f._y - op_i._y
    ndz = np_f._z - op_i._z

    dx = ndx - odx
    dy = ndy - ody
    dz = ndz - odz
    d = sqrt(dx*dx+dy*dy+dz*dz)

    err= d/org_d*100.0
    @test err==0.0
end

test_backdynamics()

function test_grad(;npart=500, nsub=50, stepc=10)
    ixs = randperm(npart)[1:nsub]

    sim = Simulation(worldnormal(npart, smth=0.01, opening_alpha=0.7),
                limit_by_steps=true, stepc=stepc)

    # getting fast gradients
    opt = Optimization(sim)
    @time exec!(sim, opt; silent=true, accumulate_history=true)
    @time grad!(opt,sim)
    fast_gx = opt.gx[ixs]

    # getting slow gradients
    slow_gx = zeros(nsub)
    exec!(sim; silent=true)
    D = 1.e-7
    g0 = grade(sim, opt)
    for i in 1:nsub
        ix = ixs[i]
        sim.xi[ix] += D
        exec!(sim; silent=true)
        sim.xi[ix] -= D
        g1 = grade(sim, opt)
        slow_gx[i] = (g1-g0)/D
    end

    fast_gx, slow_gx
end

test_grad()
