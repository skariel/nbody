immutable CompiledTreeNode
    cm_x::Float64
    cm_y::Float64
    cm_z::Float64
    m::Float64
    l::Float64
    next::Int64
end

type CompiledTree
    tree::Array{CompiledTreeNode, 1}
    nodes_used::Int64
    CompiledTree(n::Int) = new(Array(CompiledTreeNode, 2*n), 0)
end

@inline CompiledTreeNode(n::OctTreeNode{Particle}, ct::CompiledTree) =
    CompiledTreeNode(n.point._x, n.point._y, n.point._z,
        n.point._m, 2.0*n.r, isleaf(n) ? -1 : ct.nodes_used+1)

@inline function stop_cond(q::OctTreeNode{Particle}, ct::CompiledTree)
    isemptyleaf(q) && return true # empty node, nothing to do
    ct.nodes_used += 1
    @inbounds ct.tree[ct.nodes_used] = CompiledTreeNode(q, ct)
    return false
end

function compile(ct::CompiledTree, t::OctTree{Particle})
    ct.nodes_used = 0
    map(t, ct)
    ct
end
