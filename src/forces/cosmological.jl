# Post processing the newtonian stuff

@inline function calculate_accel_on_particle!(w::World{Cosmological}, i::Int64)
    __calculate_accel_on_particle!(w,i)
    @inline w.ax[i] += FAC1*w.particles[i]._x
    @inline w.ay[i] += FAC1*w.particles[i]._y
    @inline w.az[i] += FAC1*w.particles[i]._z
end

@inline function calc_accel!(w::World{Cosmological}, rng::UnitRange{Int64})
    __calc_accel!(w, rng)
    @inline for i in rng
        w.ax[i] += FAC1*w.particles[i]._x
        w.ay[i] += FAC1*w.particles[i]._y
        w.az[i] += FAC1*w.particles[i]._z
    end
end

@inline function calc_accel!(w::World{Cosmological},
    tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1},
    tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1},
    w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    __calc_accel!(w,tx,ty,tz,tax,tay,taz,w_rng,t_rng)
    @inline for i in w_rng
        w.ax[i] += FAC1*w.particles[i]._x
        w.ay[i] += FAC1*w.particles[i]._y
        w.az[i] += FAC1*w.particles[i]._z
    end
    @inline for i in t_rng
        tax[i] += FAC1*tx[i]
        tay[i] += FAC1*ty[i]
        taz[i] += FAC1*tz[i]
    end
end

@inline function calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Cosmological})
    ax, ay, az = __calc_accel(p, tree, w)
    ax += FAC1*p._x
    ay += FAC1*p._y
    az += FAC1*p._z
end
