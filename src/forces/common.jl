@inline function modify(q::OctTreeNode{Particle}, p::Particle)
    const total_mass = q.point._m + p._m
    const newx = (q.point._x*q.point._m + p._x*p._m)/total_mass
    const newy = (q.point._y*q.point._m + p._y*p._m)/total_mass
    const newz = (q.point._z*q.point._m + p._z*p._m)/total_mass
    q.point = Particle(newx, newy, newz, total_mass)
    nothing
end

@inline createtree(w::World) = OctTree(Particle; n=int(3.9*w.n))

function buildtree!(w::World, tree::OctTree{Particle})
    clear!(tree)
    # calculate new boundries same extent on both x and y
    minc = Float64(1.e30)
    maxc = Float64(-1.e30)
    for i in 1:w.n
        const p = w.particles[i]
        if p._x < minc
            minc = p._x
        end
        if p._y < minc
            minc = p._y
        end
        if p._z < minc
            minc = p._z
        end
        if p._x > maxc
            maxc = p._x
        end
        if p._y > maxc
            maxc = p._y
        end
        if p._z > maxc
            maxc = p._z
        end
    end
    r = 0.5*(maxc-minc)
    md= 0.5*(maxc+minc)
    initnode!(tree.head, r*1.05, md, md, md)
    insert!(tree, deepcopy(w.particles), Modify)
    compile!(w.tree, tree)
    nothing
end

function calc_accel!(w::World)
    chunks = get_chunks(w.n)
    @sync for i in 1:length(workers())
        @async remotecall_wait(workers()[i], calc_accel!, w, chunks[i])
    end
end

function calc_accel!(w::World, tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1}, tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1})
    w_chunks = get_chunks(w.n)
    t_chunks = get_chunks(length(tx))
    @sync for i in 1:length(workers())
        @async remotecall_wait(workers()[i], calc_accel!, w, tx,ty,tz,tax,tay,taz, w_chunks[i], t_chunks[i])
    end
end
