type Optimization
    x0::SharedArray{Float64,1}
    y0::SharedArray{Float64,1}
    z0::SharedArray{Float64,1}
    gx::SharedArray{Float64,1}
    gy::SharedArray{Float64,1}
    gz::SharedArray{Float64,1}
    sim::Simulation
    grad::SharedArray{Float64,1}
    function Optimization(sim::Simulation)
        x0 = SharedArray(Float64, sim.w.n)
        y0 = SharedArray(Float64, sim.w.n)
        z0 = SharedArray(Float64, sim.w.n)
        gx = SharedArray(Float64, sim.w.n)
        gy = SharedArray(Float64, sim.w.n)
        gz = SharedArray(Float64, sim.w.n)
        grad = SharedArray(Float64, sim.w.n)
        for i in 1:sim.w.n
             x0[i] = sim.xi[i]
             y0[i] = sim.yi[i]
             z0[i] = sim.zi[i]
             gx[i] = 0.0
             gy[i] = 0.0
             gz[i] = 0.0
             grad[i] = 0.0
        end
        new(x0,y0,z0, gx,gy,gz, sim, grad)
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
    fa = 0.2
    fb = 0.55
    a = 1.e30 # whatever
    b = 1.e30 # whatever
    recalc_a=true
    recalc_b=true
    while true
        iter += 1

        if recalc_a
            a = grade(opt, fa)
            recalc_a = false
        end
        if recalc_b
            b = grade(opt, fb)
            recalc_b = false
        end
        iter >= 12 && break

        if a > curr_grade
            fb = fa
            b = a
            fa /= 1.5
            recalc_a = true
        elseif b < a
            fa = fb
            a = b
            fb *= 1.45
            recalc_b = true
        else
            fb = fa
            b = a
            fa /= 1.4
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
        println("SHAKING!!!")
        for i in 1:opt.sim.w.n
            opt.gx[i] = randn()
            opt.gy[i] = randn()
            opt.gz[i] = randn()
        end
        return search!(opt, true)
    end
    println("   - ls: iter=",iter," dgrade=",dgrade, " f=", f)
    f
end

function optimize(opt::Optimization, maxstep=10, ming=0.001)
    step = 1
    g = 1.e30 # infinity, ha!
    f = 1.0
    fa = Float64[]
    while true
        try
            exec!(opt.sim, false, silent=true)
            g = grade(opt)
        catch e
        end

        println("step=",step," grade=",g," n=",int(f*opt.sim.w.n))
        if step >= maxstep || g <= ming
            break
        end

        # getting the gradient
        grad!(opt)

        # filtering gradient to worst 10%
        ga = sort(opt.grad)
        th = ga[int(f*length(ga))]
        for i in 1:opt.sim.w.n
            opt.grad[i] < th && continue
            opt.gx[i] = 0.0
            opt.gy[i] = 0.0
            opt.gz[i] = 0.0
        end
        f *= 0.9
        if int(f*length(ga)) < 1
            break
            f = 1.0
        end


        # linear search
        fact = 0.0001
        try
            fact = search!(opt)
            push!(fa, fact)
        catch e
        end
        if fact < 0
            break
        end

        # evolve the ics
        evolve_ics!(opt, fact)

        step += 1
    end
    g
end
