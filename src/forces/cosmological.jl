# Post processing the newtonian stuff

function calculate_accel_on_particle!(w::World{Cosmological}, i::Int64)
    __calculate_accel_on_particle!(w,i)
    F1 = FAC1(w)
    w.ax[i] += F1*w.particles[i]._x
    w.ay[i] += F1*w.particles[i]._y
    w.az[i] += F1*w.particles[i]._z
end

function calc_accel!(w::World{Cosmological}, rng::UnitRange{Int64})
    __calc_accel!(w, rng)
    F1 = FAC1(w)
    for i in rng
        w.ax[i] += F1*w.particles[i]._x
        w.ay[i] += F1*w.particles[i]._y
        w.az[i] += F1*w.particles[i]._z
    end
end

function calc_accel!(w::World{Cosmological},
    tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1},
    tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1},
    w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    __calc_accel!(w,tx,ty,tz,tax,tay,taz,w_rng,t_rng)
    F1 = FAC1(w)
    for i in w_rng
        w.ax[i] += F1*w.particles[i]._x
        w.ay[i] += F1*w.particles[i]._y
        w.az[i] += F1*w.particles[i]._z
    end
    if t_rng.start > 0
        for i in t_rng
            tax[i] += F1*tx[i]
            tay[i] += F1*ty[i]
            taz[i] += F1*tz[i]
        end
    end
end

function calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Cosmological})
    ax, ay, az = __calc_accel(p, tree, w)
    F1 = FAC1(w)
    ax += F1*p._x
    ay += F1*p._y
    az += F1*p._z
end
