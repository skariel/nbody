using Winston
import Winston: plot

plot(w::World) =
    scatter(
        [p._x for p in w.particles],
        [p._y for p in w.particles],
        0.0001, aspect_ratio=1
    )

function plot(sim::Simulation; s1=0.01, s2=0.02, color1="red", color2="black")
    hold(false)
    scatter(
        [p._x for p in sim.w.particles],
        [p._y for p in sim.w.particles],
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
