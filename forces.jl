@inline function modify(q::OctTreeNode{Particle}, p::Particle)
    const total_mass = q.point._m + p._m
    const newx = (q.point._x*q.point._m + p._x*p._m)/total_mass
    const newy = (q.point._y*q.point._m + p._y*p._m)/total_mass
    const newz = (q.point._z*q.point._m + p._z*p._m)/total_mass
    q.point = Particle(newx, newy, newz, total_mass)
end

function buildtree(w::World)
    clear!(w.tree)
    # calculate new boundries same extent on both x and y
    minc = Float64(1.e30)
    maxc = Float64(-1.e30)
    for i in 1:w.n
        @inbounds const p = w.particles[i]
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
    initnode!(w.tree.head, r*1.05, md, md, md)
    insert!(w.tree, deepcopy(w.particles), Modify)
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

@inline function stop_cond(q::OctTreeNode{Particle}, data::DataToCalculateAccelOnParticle)
    isemptyleaf(q) && return true # empty node, nothing to do

    if isleaf(q)
        # we have a single particle in the node
        q.point._x == data.px &&
        q.point._y == data.py &&
        q.point._z == data.pz && return true
        const dx = q.point._x - data.px
        const dx2 = dx*dx
        const dy = q.point._y - data.py
        const dy2 = dy*dy
        const dz = q.point._z - data.pz
        const dz2 = dz*dz
        const dr2 = dx2+dy2+dz2+data.w.smth2
        const dr = sqrt(dr2)
        const denom = dr2*dr/q.point._m
        data.ax += dx/denom
        data.ay += dy/denom
        data.az += dz/denom
        return true
    end

    # here q is divided. Check we are not too close to the cell
    const lx = 2.0*q.r
    const lx2 = lx*lx
    const dqx = q.midx - data.px
    const dqx2 = dqx*dqx
    const dqy = q.midy - data.py
    const dqy2 = dqy*dqy
    const dqz = q.midz - data.pz
    const dqz2 = dqz*dqz
    const fac = data.w.opening_excluded_frac2*lx2
    dqx2 < fac && dqy2 < fac && dqz2 < fac && return false # we need to further open the node
    const dx = q.point._x - data.px
    const dx2 = dx*dx
    const dy = q.point._y - data.py
    const dy2 = dy*dy
    const dz = q.point._z - data.pz
    const dz2 = dz*dz
    const r2 = dx2 + dy2 + dz2
    lx2/r2 > data.w.opening_alpha2 && return false # we need to further open the node

    # consider the node, no further need to open it
    const dr2 = r2+data.w.smth2
    const dr = sqrt(dr2)
    const denom = dr2*dr/q.point._m
    data.ax += dx/denom
    data.ay += dy/denom
    data.az += dz/denom
    return true
end

@inline function calc_accel(p::Particle, tree::OctTree{Particle}, w::World)
    @inbounds data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,p._x,p._y,p._z,w)
    map(tree, data)
    data.ax, data.ay, data.az
end

@inline function calculate_accel_on_particle!(w::World, particle_ix::Int64)
    @inbounds const p = w.particles[particle_ix]
    @inbounds data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,p._x,p._y,p._z,w)
    map(w.tree, data)
    @inbounds w.ax[particle_ix] = data.ax
    @inbounds w.ay[particle_ix] = data.ay
    @inbounds w.az[particle_ix] = data.az
end

function calc_accel!(w::World)
    buildtree(w)
    data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in 1:w.n
        const p = w.particles[i]
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
    end
end

#################################################################
#
#  Forces on compiled trees
#
#################################################################

@inline function stop_cond(q::CompiledTreeNode, data::DataToCalculateAccelOnParticle)
    if q.l < 0
        # this is a leaf
        # we have a single particle in the node
        q.cm_x == data.px &&
        q.cm_y == data.py &&
        q.cm_z == data.pz && return true
        const dx = q.cm_x - data.px
        const dx2 = dx*dx
        const dy = q.cm_y - data.py
        const dy2 = dy*dy
        const dz = q.cm_z - data.pz
        const dz2 = dz*dz
        const dr2 = dx2+dy2+dz2+data.w.smth2
        const dr = sqrt(dr2)
        const denom = dr2*dr/q.m
        data.ax += dx/denom
        data.ay += dy/denom
        data.az += dz/denom
        return true
    end

    # here q is divided. Check we are not too close to the cell
    const l2 = q.l*q.l
    const dx = q.cm_x - data.px
    const dx2 = dx*dx
    const dy = q.cm_y - data.py
    const dy2 = dy*dy
    const dz = q.cm_z - data.pz
    const dz2 = dz*dz
    const r2 = dx2 + dy2 + dz2
    l2/r2 > data.w.opening_alpha2 && return false # we need to further open the node

    # consider the node, no further need to open it
    const dr2 = r2+data.w.smth2
    const dr = sqrt(dr2)
    const denom = dr2*dr/q.m
    data.ax += dx/denom
    data.ay += dy/denom
    data.az += dz/denom
    return true
end

function calc_accel!(w::World, ::Type{Compile})
    ct = CompiledTree(w.n)
    calc_accel!(w, Compile, ct)
end

function calc_accel!(w::World, ::Type{Compile}, ct::CompiledTree)
    buildtree(w)
    compile(ct, w.tree)
    data = DataToCalculateAccelOnParticle(0.0,0.0,0.0,0.0,0.0,0.0,w)
    @inbounds for i in 1:w.n
        const p = w.particles[i]
        data.ax = 0.0
        data.ay = 0.0
        data.az = 0.0
        data.px = p._x
        data.py = p._y
        data.pz = p._z
        map(ct, data)
        @inbounds w.ax[i] = data.ax
        @inbounds w.ay[i] = data.ay
        @inbounds w.az[i] = data.az
    end
end
