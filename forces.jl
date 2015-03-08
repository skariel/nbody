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

type DataToCalculateNewtonianAccelOnParticle
    ax::Float64 # to store calculated acceleration
    ay::Float64 # ...
    az::Float64 # ...
    px::Float64 # particle position
    py::Float64 # ...
    pz::Float64 # ...
    w::World{Newtonian}
end

type DataToCalculateCosmologicalAccelOnParticle
    ax::Float64 # to store calculated acceleration
    ay::Float64 # ...
    az::Float64 # ...
    px::Float64 # particle position
    py::Float64 # ...
    pz::Float64 # ...
    vx::Float64 # particle velocity
    vy::Float64 # ..
    vz::Float64 # ..
    w::World{Cosmological}
end

@inline function stop_cond(q::CompiledOctTreeNode{Particle}, data::DataToCalculateNewtonianAccelOnParticle)
    const dx = q.point._x - data.px
    const dy = q.point._y - data.py
    const dz = q.point._z - data.pz

    q.l < 0.0 && dx==0.0 && dy==0.0 && dz==0.0 && return true

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    dr2 < 1.e-6*1.e-6 && return true

    q.l > 0 && q.l*q.l/dr2 > data.w.opening_alpha2 && return false # we need to further open the node

    const smthdr2 = dr2+data.w.smth2
    const smthdr = sqrt(smthdr2)
    const denom = smthdr2*smthdr/q.point._m
    data.ax += dx/denom
    data.ay += dy/denom
    data.az += dz/denom
    return true
end

@inline function stop_cond(q::CompiledOctTreeNode{Particle}, data::DataToCalculateCosmologicalAccelOnParticle)
    const dx = q.point._x - data.px
    const dy = q.point._y - data.py
    const dz = q.point._z - data.pz

    q.l < 0.0 && dx==0.0 && dy==0.0 && dz==0.0 && return true

    const dx2 = dx*dx
    const dy2 = dy*dy
    const dz2 = dz*dz
    const dr2 = dx2 + dy2 + dz2

    dr2 < 1.e-6*1.e-6 && return true

    q.l > 0 && q.l*q.l/dr2 > data.w.opening_alpha2 && return false # we need to further open the node

    const smthdr2 = dr2+data.w.smth2
    const smthdr = sqrt(smthdr2)
    const denom = smthdr2*smthdr/q.point._m*data.w.space.a3
    data.ax += dx/denom - data.w.space.adda*dx - data.w.space.H2*data.vx
    data.ay += dy/denom - data.w.space.adda*dy - data.w.space.H2*data.vy
    data.az += dz/denom - data.w.space.adda*dz - data.w.space.H2*data.vz
    return true
end

function calculate_accel_on_particle!(w::World{Newtonian}, i::Int64)
    data = DataToCalculateNewtonianAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds const p = w.particles[i]
    data.ax = 0.0
    data.ay = 0.0
    data.az = 0.0
    data.px = p._x
    data.py = p._y
    data.pz = p._z
    map(w.tree, data)
    @inbounds w.ax[i] = data.ax
    @inbounds w.ay[i] = data.ay
    @inbounds w.az[i] = data.az
    nothing
end

function calc_accel!(w::World{Newtonian}, rng::UnitRange{Int64})
    data = DataToCalculateNewtonianAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
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

function calc_accel!(w::World{Newtonian}, tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1}, tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1}, w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    data = DataToCalculateNewtonianAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in w_rng
        p = w.particles[i]
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
    @inbounds for i in t_rng
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.px = tx[i]
        data.py = ty[i]
        data.pz = tz[i]
        map(w.tree, data)
        tax[i] = data.ax
        tay[i] = data.ay
        taz[i] = data.az
    end
    nothing
end

function calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Newtonian})
    data = DataToCalculateNewtonianAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    data.ax = 0.0
    data.ay = 0.0
    data.az = 0.0
    data.px = p._x
    data.py = p._y
    data.pz = p._z
    map(tree, data)
    data.ax, data.ay, data.az
end

### Cosmological forces

function calculate_accel_on_particle!(w::World{Cosmological}, i::Int64)
    data = DataToCalculateCosmologicalAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,w)
    const p = w.particles[i]
    data.ax = 0.0
    data.ay = 0.0
    data.az = 0.0
    @inbounds data.vx = w.vx[i]
    @inbounds data.vy = w.vy[i]
    @inbounds data.vz = w.vz[i]
    data.px = p._x
    data.py = p._y
    data.pz = p._z
    map(w.tree, data)
    @inbounds w.ax[i] = data.ax
    @inbounds w.ay[i] = data.ay
    @inbounds w.az[i] = data.az
    nothing
end

function calc_accel!(w::World{Cosmological}, rng::UnitRange{Int64})
    data = DataToCalculateCosmologicalAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in rng
        const p = w.particles[i]
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.vx = w.vx[i]
        data.vy = w.vy[i]
        data.vz = w.vz[i]
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

function calc_accel!(w::World{Cosmological}, tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1}, tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1}, w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    data = DataToCalculateCosmologicalAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in w_rng
        p = w.particles[i]
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.vx = w.vx[i]
        data.vy = w.vy[i]
        data.vz = w.vz[i]
        data.px = p._x
        data.py = p._y
        data.pz = p._z
        map(w.tree, data)
        w.ax[i] = data.ax
        w.ay[i] = data.ay
        w.az[i] = data.az
    end
    @inbounds for i in t_rng
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.vx = w.vx[i]
        data.vy = w.vy[i]
        data.vz = w.vz[i]
        data.px = tx[i]
        data.py = ty[i]
        data.pz = tz[i]
        map(w.tree, data)
        tax[i] = data.ax
        tay[i] = data.ay
        taz[i] = data.az
    end
    nothing
end

function calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Cosmological})
    data = DataToCalculateCosmologicalAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,w)
    data.ax = 0.0
    data.ay = 0.0
    data.az = 0.0
    @inbounds data.vx = w.vx[i]
    @inbounds data.vy = w.vy[i]
    @inbounds data.vz = w.vz[i]
    data.px = p._x
    data.py = p._y
    data.pz = p._z
    map(tree, data)
    data.ax, data.ay, data.az
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
