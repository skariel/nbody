# Post processing the newtonian stuff

@inline function calculate_accel_on_particle!(w::World{Cosmological}, i::Int64)
    __calculate_accel_on_particle!(w,i)
    @inline w.ax[i] /= w.space.a
    @inline w.ay[i] /= w.space.a
    @inline w.az[i] /= w.space.a
end

@inline function calc_accel!(w::World{Cosmological}, rng::UnitRange{Int64})
    __calc_accel!(w, rng)
    @inline for i in rng
        w.ax[i] /= w.space.a
        w.ay[i] /= w.space.a
        w.az[i] /= w.space.a
    end
end

@inline function calc_accel!(w::World{Cosmological},
    tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1},
    tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1},
    w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    __calc_accel!(w,tx,ty,tz,tax,tay,taz,w_rng,t_rng)
    @inline for i in w_rng
        w.ax[i] /= w.space.a
        w.ay[i] /= w.space.a
        w.az[i] /= w.space.a
    end
    @inline for i in t_rng
        tax[i] /= w.space.a
        tay[i] /= w.space.a
        taz[i] /= w.space.a
    end
end

@inline function calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Cosmological})
    ax, ay, az = __calc_accel(p, tree, w)
    ax /= w.space.a
    ay /= w.space.a
    az /= w.space.a
end
