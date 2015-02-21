@inline function stop_cond(q::OctTreeNode{Particle}, ct::CompiledTree)
    isemptyleaf(q) && return true # empty node, nothing to do
    ct.nodes_used += 1
    @inbounds ct.tree[ct.nodes_used] = CompiledTreeNode(q, ct)
    q.id = ct.nodes_used
    return false
end

function compile(ct::CompiledTree, t::OctTree{Particle})
    ct.nodes_used = 0
    map(t, ct)

    # fix neighbours
    for i in 1:t.number_of_nodes_used
        @inbounds const q=t.nodes[i]
        !q.is_divided && continue
        childs = (q.lxlylz, q.lxlyhz, q.lxhylz, q.lxhyhz, q.hxlylz, q.hxlyhz, q.hxhylz, q.hxhyhz)
        for a in 1:7
            @inbounds const qa = childs[a]
            qa.id <= 0 && continue
            # searching qa neighbor...
            qb = qa
            for b in (a+1):8
                @inbounds if childs[b].id > 0
                    qb = childs[b]
                    break
                end
            end
            is(qb,qa) && continue # no neighbor found lets leave next=-1
            # qa and qb are neighbors -- fix that in the compiled tree
            @inbounds ct.tree[qb.id] = withnext(ct.tree[qb.id], qa.id)
        end
    end
    ct
end

@inline function map(t::CompiledTree, cond_data)
    curr_stack_ix = 1
    @inbounds t.faststack[1] = 1
    while curr_stack_ix > 0
        @inbounds node_ix = t.faststack[curr_stack_ix]
        curr_stack_ix -= 1
        @inbounds if node_ix>0 && node_ix<=t.nodes_used && !stop_cond(t.tree[node_ix], cond_data)
            curr_stack_ix += 1
            next_ix = node_ix+1
            @inbounds t.faststack[curr_stack_ix] = next_ix
            while next_ix > 0
                curr_stack_ix += 1
                @inbounds next_ix = t.tree[next_ix].next
                @inbounds t.faststack[curr_stack_ix] = next_ix
            end
        end
    end
end
