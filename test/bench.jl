@everywhere begin
    import nbody
    reload("nbody")
    nb = nbody;
end

w = nb.worldnormal(1000000, smth=0.0, opening_alpha=0.7);
@show typeof(w)

tree = nb.createtree(w);

@time nb.buildtree!(w, tree)

@show length(tree.nodes)
@show tree.number_of_nodes_used
@show w.tree.number_of_nodes_used

@time nb.compile!(w.tree, tree);

@time nb.calc_accel!(w)
