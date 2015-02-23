@inline function modify(q::OctTreeNode{Particle}, p::Particle)
    const total_mass = q.point._m + p._m
    const newx = (q.point._x*q.point._m + p._x*p._m)/total_mass
    const newy = (q.point._y*q.point._m + p._y*p._m)/total_mass
    const newz = (q.point._z*q.point._m + p._z*p._m)/total_mass
    q.point = Particle(newx, newy, newz, total_mass)
    nothing
end

@inline createtree(w::World) = OctTree(Particle; n=w.n)

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

type DataToCalculateAccelOnParticle
    ax::Float64
    ay::Float64
    az::Float64
    px::Float64
    py::Float64
    pz::Float64
    w::World
end

@inline function stop_cond(q::CompiledOctTreeNode{Particle}, data::DataToCalculateAccelOnParticle)
    const dx = q.point._x - data.px
    const dy = q.point._y - data.py
    const dz = q.point._z - data.pz

    q.l < 0.0 && dx==0.0 && dy==0.0 && dz==0.0 && return true

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    q.l > 0 && q.l*q.l/dr2 > data.w.opening_alpha2 && return false # we need to further open the node

    const smthdr2 = dr2+data.w.smth2
    const smthdr = sqrt(smthdr2)
    const denom = smthdr2*smthdr/q.point._m
    data.ax += dx/denom
    data.ay += dy/denom
    data.az += dz/denom
    return true
end

function calculate_accel_on_particle!(w::World, i::Int64)
    data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    const p = w.particles[i]
    data.ax = 0.0
    data.ay = 0.0
    data.az = 0.0
    data.px = p._x
    data.py = p._y
    data.pz = p._z
    map(w.tree, data)
    w.ax[i] = data.ax
    w.ay[i] = data.ay
    w.az[i] = data.az
    nothing
end

function calc_accel!(w::World, rng::UnitRange{Int64})
    @show rng
    data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in rng
        const p = w.particles[i]
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.px = p._x
        data.py = p._y
        data.pz = p._z
        map(w.tree, data)
        w.ax[i] = data.ax
        w.ay[i] = data.ay
        w.az[i] = data.az
    end
    nothing
end

function get_chunks(n::Int64)
    chunks_i = int(linspace(1,n,1+length(workers())))
    chunks_f = chunks_i[2:end]
    chunks_f[1:end-1] -= 1
    chunks_i = chunks_i[1:end-1]
    [chunks_i[i]:chunks_f[i] for i in 1:length(chunks_i)]
end

function calc_accel!(w::World)
    chunks = get_chunks(w.n)
    @sync for i in 1:length(workers())
        @async remotecall_wait(workers()[i], calc_accel!, w, chunks[i])
    end
end
