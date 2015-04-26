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
    fa = 0.07
    fb = 0.44
    fa = 0.01
    fb = 0.07
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
        iter >= 8 && break

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


function fix_away_particles!(opt::Optimization)
    for i in 1:opt.sim.w.n
        opt.grad[i] < 0.99*opt.sim.w.smth2 && continue

        # search for closest particle0
        ixd2 = 1.e30
        ix = -1
        for j in 1:opt.sim.w.n
            dx = opt.sim.w.particles[j]._x - opt.x0[i]
            dy = opt.sim.w.particles[j]._y - opt.y0[i]
            dz = opt.sim.w.particles[j]._z - opt.z0[i]
            d2 = dx*dx+dy*dy+dz*dz
            if d2 < ixd2
                ix = j
                ixd2 = d2
            end
        end

        # se new location
        opt.sim.xi[i] = opt.sim.xi[ix] + randn()*sqrt(opt.sim.w.smth2)/5.
        opt.sim.yi[i] = opt.sim.yi[ix] + randn()*sqrt(opt.sim.w.smth2)/5.
        opt.sim.zi[i] = opt.sim.zi[ix] + randn()*sqrt(opt.sim.w.smth2)/5.
    end
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

        fix_away_particles!(opt)

        println("step=",step," grade=",g," n=",int(f*opt.sim.w.n))
        if step >= maxstep || g <= ming
            break
        end

        try
            # getting the gradient
            grad!(opt)
        catch e
        end

        # filtering gradient random 5%
        for i in 1:opt.sim.w.n
            rand() > (1.0-f) && continue
            opt.gx[i] = 0.0
            opt.gy[i] = 0.0
            opt.gz[i] = 0.0
        end
        f = sqrt(rand())
        if f < 0.07
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
