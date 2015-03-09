
type DataToCalculateNewtonianAccelOnParticle
    ax::Float64 # to store calculated acceleration
    ay::Float64 # ...
    az::Float64 # ...
    px::Float64 # particle position
    py::Float64 # ...
    pz::Float64 # ...
    w::World{Newtonian}
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
