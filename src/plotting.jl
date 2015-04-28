using Winston
import Winston: plot

function plot(w::World; sz=0.0001, color="black")
    s = scatter(
        [p._x/1000.0 for p in w.particles],
        [p._y/1000.0 for p in w.particles],
        sz, aspect_ratio=1, color=color
    )
    xlabel("Mpc/h")
    ylabel("Mpc/h")
    s
end

function plotz(w::World; sz=0.0001, color="black")
    to_z_space!(w)
    plot(w, sz=sz, color=color)
    to_c_space!(w)
    title("Redshift space")
end

function plot(sim::Simulation; s1=0.01, s2=0.02, color1="red", color2="black")
    hold(false)
    plot(sim.w, sz=s1, color=color1)
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
    plot(w, sz=s1, color=color1)
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
    plot(w1, sz=s1, color=color1)
    hold(true)
    s = plot(w2, sz=s2, color=color2)
    hold(false)
    s
end

function plot(opt::Optimization)
    hold(false)
    plot(w, sz=0.03, color="red")
    hold(true)
    scatter(
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
