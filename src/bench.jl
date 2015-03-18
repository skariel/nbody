push!(LOAD_PATH,"./src/")
@everywhere begin
    import nbody
    reload("nbody")
    nb = nbody;
end

N = 2048
stepc = 120
w = nb.worldspherical(N, smth=0.2, opening_alpha=0.5, space=nb.Cosmological, scale=170.)
sim = nb.Simulation(w, limit_by_steps=false, stepc=36, ti=0.02, tf=0.04, n_test_particle=0)
