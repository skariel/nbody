function evolve_ics!(opt::Optimization, sim::Simulation, fact::Float64)
    @inbounds for i in 1:sim.w.n
        sim.xi[i] -= opt.gx[i]*fact
        sim.yi[i] -= opt.gy[i]*fact
        sim.zi[i] -= opt.gz[i]*fact
    end
end

function grade(opt::Optimization, sim::Simulation, fact::Float64)
    evolve_ics!(opt, sim, fact)
    exec!(sim, silent=true)
    res = grade(opt, sim)
    evolve_ics!(opt, sim, -fact)
    res
end

function search!(opt::Optimization, sim::Simulation, shaking=false)
    curr_grade = grade(opt, sim)
    iter = 0
    fa = 0.3
    fb = 0.7
    a = 1.e30 # whatever
    b = 1.e30 # whatever
    recalc_a=true
    recalc_b=true
    while iter < 8
        iter += 1

        if recalc_a
            a = grade(opt, sim, fa)
            recalc_a = false
        end
        if recalc_b
            b = grade(opt, sim, fb)
            recalc_b = false
        end

        if a > curr_grade
            fb = fa
            b = a
            fa /= 1.5
            recalc_a = true
        elseif b < a
            fa = fb
            a = b
            fb *= 1.3
            recalc_b = true
        else
            fb = fa
            b = a
            fa /= 1.2
            recalc_a = true
        end


    end
    if recalc_a
        a = grade(opt, sim, fa)
    end
    if recalc_b
        b = grade(opt, sim, fb)
    end
    (g,f) = a<b ? (a,fa):(b,fb)
    dgrade = g/curr_grade
    if !shaking && dgrade > 0.9
        # shake!
        println("Shaking!")
        @inbounds for i in 1:sim.w.n
            opt.gx[i] = randn()
            opt.gy[i] = randn()
            opt.gz[i] = randn()
        end
        return search!(opt,sim,true)
    end
    println("   - ls: iter=",iter," dgrade=",dgrade, " f=", f)
    f
end

function optimize(opt::Optimization, sim::Simulation, maxstep=10, ming=0.001)
    step = 1
    g = 1.e30 # infinity, ha!
    while step <= maxstep && g > ming
        exec!(sim, silent=true)
        g = grade(opt,sim)
        println("step=",step," grade=",g)

        # getting the gradient
        grad!(opt, sim)

        # linear search
        fact = search!(opt, sim)

        # evolve the ics
        evolve_ics!(opt, sim, fact)

        step += 1
    end
end
