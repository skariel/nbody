using Winston
import Winston: plot

plot(w::World) =
    scatter(
        [p._x for p in w.particles],
        [p._y for p in w.particles],
        0.00001
    )
