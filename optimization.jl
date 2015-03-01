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

function grade_rnd(opt::Optimization, sim::Simulation, fact::Float64)
    rx = randn(sim.w.n)*fact
    ry = randn(sim.w.n)*fact
    rz = randn(sim.w.n)*fact
    sim.xi += rx
    sim.yi += ry
    sim.zi += rz
    exec!(sim, silent=true)
    res = grade(opt, sim)
    sim.xi -= rx*fact
    sim.yi -= ry*fact
    sim.zi -= rz*fact
    res
end

function search!(opt::Optimization, sim::Simulation)
    curr_grade = grade(opt, sim)
    iter = 0
    a = 1.0e30
    b = 1.0e30
    fa = 0.0
    fb = 0.1
    while a>curr_grade && iter < 100
        fa = 0.0
        fb = 0.1
        while fb-fa > 1.e-4
            iter += 1
            a = grade(opt, sim, fa)
            b = grade(opt, sim, fb)
            if b < a
                fa = fb
                fb *= 2.0
            else
                fb = 0.5*(fb+fa)
            end
        end
        if a > curr_grade
            a = grade_rnd(opt, sim, std(sim.xi)/500.0)
        end
    end
    println("   - ls: iter=",iter," dgrade=",a/curr_grade, " fa=", fa)
    fa
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
