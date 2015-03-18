type Optimization
    x0::SharedArray{Float64,1}
    y0::SharedArray{Float64,1}
    z0::SharedArray{Float64,1}
    gx::SharedArray{Float64,1}
    gy::SharedArray{Float64,1}
    gz::SharedArray{Float64,1}
    sim::Simulation
    function Optimization(sim::Simulation)
        x0 = SharedArray(Float64, sim.w.n)
        y0 = SharedArray(Float64, sim.w.n)
        z0 = SharedArray(Float64, sim.w.n)
        gx = SharedArray(Float64, sim.w.n)
        gy = SharedArray(Float64, sim.w.n)
        gz = SharedArray(Float64, sim.w.n)
        for i in 1:sim.w.n
             x0[i] = sim.xi[i]
             y0[i] = sim.yi[i]
             z0[i] = sim.zi[i]
             gx[i] = 0.0
             gy[i] = 0.0
             gz[i] = 0.0
        end
        new(x0,y0,z0, gx,gy,gz, sim)
    end
end

function evolve_ics!(opt::Optimization, fact::Float64)
     for i in 1:opt.sim.w.n
        opt.sim.xi[i] -= opt.gx[i]*fact
        opt.sim.yi[i] -= opt.gy[i]*fact
        opt.sim.zi[i] -= opt.gz[i]*fact
    end
end

function grade(opt::Optimization, fact::Float64)
    evolve_ics!(opt, fact)
    exec!(opt.sim, false, silent=true)
    res = grade(opt)
    evolve_ics!(opt, -fact)
    res
end

function search!(opt::Optimization, shaking=false)
    curr_grade = grade(opt)
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
            a = grade(opt, fa)
            recalc_a = false
        end
        if recalc_b
            b = grade(opt, fb)
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
        a = grade(opt, fa)
    end
    if recalc_b
        b = grade(opt, fb)
    end
    (g,f) = a<b ? (a,fa):(b,fb)
    dgrade = g/curr_grade
    if !shaking && dgrade > 0.97
        # shake!
        println("Shaking!")
         for i in 1:opt.sim.w.n
            opt.gx[i] = randn()
            opt.gy[i] = randn()
            opt.gz[i] = randn()
        end
        return search!(opt,true)
    end
    println("   - ls: iter=",iter," dgrade=",dgrade, " f=", f)
    f
end

function optimize(opt::Optimization, maxstep=10, ming=0.001)
    step = 1
    g = 1.e30 # infinity, ha!
    while step <= maxstep && g > ming
        exec!(opt.sim, false, silent=true)
        g = grade(opt)
        println("step=",step," grade=",g)

        # getting the gradient
        grad!(opt)

        # linear search
        fact = search!(opt)

        # evolve the ics
        evolve_ics!(opt, fact)

        step += 1
    end
end
