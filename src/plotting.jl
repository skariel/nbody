using Winston
import Winston: plot

function plot(w::World; sz=0.0001)
    scatter(
        [p._x/1000.0 for p in w.particles],
        [p._y/1000.0 for p in w.particles],
        sz, aspect_ratio=1
    )
    xlabel("Mpc/h")
    ylabel("Mpc/h")
end

function plot(sim::Simulation; s1=0.01, s2=0.02, color1="red", color2="black")
    hold(false)
    scatter(
        [p._x/1000.0 for p in sim.w.particles],
        [p._y/1000.0 for p in sim.w.particles],
        s1, aspect_ratio=1, color=color1
    )
    hold(true)
    s = scatter(
        sim.xi./1000.0, sim.yi./1000.0,
        s2, color=color2, aspect_ratio=1
    )
    xlabel("Mpc/h")
    ylabel("Mpc/h")
    hold(false)
    s
end

function plot(sim::Simulation, w::World; s1=0.01, s2=0.02, color1="red", color2="black")
    hold(false)
    scatter(
        [p._x for p in w.particles],
        [p._y for p in w.particles],
        s1, aspect_ratio=1, color=color1
    )
    hold(true)
    s = scatter(
        sim.xi, sim.yi,
        s2, color=color2, aspect_ratio=1
    )
    hold(false)
    s
end

function plot(w1::World, w2::World; s1=0.01, s2=0.02, color1="red", color2="black")
    hold(false)
    scatter(
        [p._x for p in w1.particles],
        [p._y for p in w1.particles],
        s1, aspect_ratio=1, color=color1
    )
    hold(true)
    s = scatter(
        [p._x for p in w2.particles],
        [p._y for p in w2.particles],
        s2, aspect_ratio=1, color=color2
    )
    hold(false)
    s
end

function plot(opt::Optimization)
    hold(false)
    scatter(
        [p._x/1000.0 for p in opt.sim.w.particles],
        [p._y/1000.0 for p in opt.sim.w.particles],
        0.03, aspect_ratio=1, color="red"
    )
    hold(true)
    s = scatter(
        opt.x0./1000.0,
        opt.y0./1000.0,
        0.02, aspect_ratio=1, color="black"
    )
    s = scatter(
        opt.sim.xi./1000.0,
        opt.sim.yi./1000.0,
        0.01, aspect_ratio=1, color="green"
    )
    xlabel("Mpc/h")
    ylabel("Mpc/h")
    hold(false)
    s
end
