using Base.Test

function calc_accel_brute_force(w::World, ixs=1:w.n)
    ax = zeros(w.n)
    ay = zeros(w.n)
    az = zeros(w.n)
    for i in ixs
        @inbounds p_i = w.particles[i]
        for j in 1:w.n
            i==j && continue
            @inbounds pj = w.particles[j]
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
            @inbounds ax[i] += dx/denom
            @inbounds ay[i] += dy/denom
            @inbounds az[i] += dz/denom
        end
    end
    ax[ixs],ay[ixs],az[ixs]
end

type TestAcc
    ax_tree::Array{Float64,1}
    a_tree::Array{Float64,1}
    ax_bf::Array{Float64,1}
    a_bf::Array{Float64,1}
    ferr::Array{Float64,1}
    fbelow::Array{Float64,1}
end

function test_acc(n, nout)
    ixs = randperm(n)[1:nout]
    w = worldspherical(n)
    buildtree(w)

    for i in ixs
        calculate_accel_on_particle!(w, i)
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

function test()
    w = worldnormal(10000)
    buildtree(w)

    #testing all points are insize tree boundaries:
    for p in w.particles
        @test p._x > w.tree.head.midx - w.tree.head.r
        @test p._x < w.tree.head.midx + w.tree.head.r
        @test p._y > w.tree.head.midy - w.tree.head.r
        @test p._y < w.tree.head.midy + w.tree.head.r
        @test p._z > w.tree.head.midz - w.tree.head.r
        @test p._z < w.tree.head.midz + w.tree.head.r
    end

    # testing number of full leafs
    tot_not_empty = 0
    tot_massive_leafs = 0
    for n in w.tree.nodes
        if isfullleaf(n)
            tot_not_empty += 1
        end
        if n.point._m > 1.e-13 && isleaf(n)
            tot_massive_leafs += 1
        end
    end
    @test tot_not_empty == w.n
    @test tot_massive_leafs == w.n

    total_mass = 0.0
    for n in w.tree.nodes[1:w.tree.number_of_nodes_used]
        if !n.is_empty
            total_mass += n.point._m
        end
    end
    @test_approx_eq_eps total_mass 1.0 1.e-4

    # testing tree construction
    for n in w.tree.nodes

        # testing children radius is 0.5*r
        if n.is_divided
            @test_approx_eq n.r/2 n.lxlylz.r
            @test_approx_eq n.r/2 n.lxlyhz.r
            @test_approx_eq n.r/2 n.lxhylz.r
            @test_approx_eq n.r/2 n.lxhyhz.r
            @test_approx_eq n.r/2 n.hxlylz.r
            @test_approx_eq n.r/2 n.hxlyhz.r
            @test_approx_eq n.r/2 n.hxhylz.r
            @test_approx_eq n.r/2 n.hxhyhz.r
        end

        # testing mass in children is mass in parent
        if n.is_divided
            parent_mass = n.point._m
            children_mass =
                n.lxlylz.point._m  +
                n.lxlyhz.point._m  +
                n.lxhylz.point._m  +
                n.lxhyhz.point._m  +
                n.hxlylz.point._m  +
                n.hxlyhz.point._m  +
                n.hxhylz.point._m  +
                n.hxhyhz.point._m
            @test_approx_eq_eps parent_mass children_mass 1.e-4
        end

        # testing all divided nodes are empty
        if n.is_divided
            @test n.is_empty
            @test !isleaf(n)
            @test !isemptyleaf(n)
            @test !isfullleaf(n)
        else
            @test isleaf(n)
            if isfullleaf(n)
                @test isfullleaf(n)
            end
        end

        # testing center of mass in children builds com in parent
        if n.is_divided
            parent_x = n.point._x
            parent_y = n.point._y
            parent_z = n.point._z
            children_x = (
                n.lxlylz.point._x*n.lxlylz.point._m  +
                n.lxlyhz.point._x*n.lxlyhz.point._m  +
                n.lxhylz.point._x*n.lxhylz.point._m  +
                n.lxhyhz.point._x*n.lxhyhz.point._m  +
                n.hxlylz.point._x*n.hxlylz.point._m  +
                n.hxlyhz.point._x*n.hxlyhz.point._m  +
                n.hxhylz.point._x*n.hxhylz.point._m  +
                n.hxhyhz.point._x*n.hxhyhz.point._m
                )/n.point._m
            children_y = (
                n.lxlylz.point._y*n.lxlylz.point._m  +
                n.lxlyhz.point._y*n.lxlyhz.point._m  +
                n.lxhylz.point._y*n.lxhylz.point._m  +
                n.lxhyhz.point._y*n.lxhyhz.point._m  +
                n.hxlylz.point._y*n.hxlylz.point._m  +
                n.hxlyhz.point._y*n.hxlyhz.point._m  +
                n.hxhylz.point._y*n.hxhylz.point._m  +
                n.hxhyhz.point._y*n.hxhyhz.point._m
                )/n.point._m
            children_z = (
                n.lxlylz.point._z*n.lxlylz.point._m  +
                n.lxlyhz.point._z*n.lxlyhz.point._m  +
                n.lxhylz.point._z*n.lxhylz.point._m  +
                n.lxhyhz.point._z*n.lxhyhz.point._m  +
                n.hxlylz.point._z*n.hxlylz.point._m  +
                n.hxlyhz.point._z*n.hxlyhz.point._m  +
                n.hxhylz.point._z*n.hxhylz.point._m  +
                n.hxhyhz.point._z*n.hxhyhz.point._m
                )/n.point._m
            @test_approx_eq_eps parent_x children_x 1.e-4
            @test_approx_eq_eps parent_y children_y 1.e-4
            @test_approx_eq_eps parent_z children_z 1.e-4

        end

    end

    println("*** All tests passed! ***")
end

test()
