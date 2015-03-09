# Post processing the newtonian stuff

@inline post(w::World{Cosmological}, x::Float64, v::Float64, a::Float64) =
    a/w.space.a3 - w.space.ddaa3*x - w.space.FU*v

@inline function calculate_accel_on_particle!(w::World{Cosmological}, i::Int64)
    __calculate_accel_on_particle!(w,i)
    @inline w.ax[i] = post(w, w.particles[i]._x, w.vx[i], w.ax[i])
    @inline w.ay[i] = post(w, w.particles[i]._y, w.vy[i], w.ay[i])
    @inline w.az[i] = post(w, w.particles[i]._z, w.vz[i], w.az[i])
end

@inline function calc_accel!(w::World{Cosmological}, rng::UnitRange{Int64})
    __calc_accel!(w, rng)
    @inline for i in rng
        w.ax[i] = post(w, w.particles[i]._x, w.vx[i], w.ax[i])
        w.ay[i] = post(w, w.particles[i]._y, w.vy[i], w.ay[i])
        w.az[i] = post(w, w.particles[i]._z, w.vz[i], w.az[i])
    end
end

@inline function calc_accel!(w::World{Cosmological},
    tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1},
    tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1},
    tvx::SharedArray{Float64, 1}, tvy::SharedArray{Float64, 1}, tvz::SharedArray{Float64, 1},
    w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64})
    __calc_accel!(w,tx,ty,tz,tax,tay,taz,w_rng,t_rng)
    @inline for i in w_rng
        w.ax[i] = post(w, w.particles[i]._x, w.vx[i], w.ax[i])
        w.ay[i] = post(w, w.particles[i]._y, w.vy[i], w.ay[i])
        w.az[i] = post(w, w.particles[i]._z, w.vz[i], w.az[i])
    end
    @inline for i in t_rng
        tax[i] = post(w, tx[i], w.vx[i], tax[i])
        tay[i] = post(w, ty[i], w.vy[i], tay[i])
        taz[i] = post(w, tz[i], w.vz[i], taz[i])
    end
end

@inline function calc_accel(p::Particle, vx::Float64, vy::Float64, vz::Float64, tree::CompiledOctTree{Particle}, w::World{Cosmological})
    ax, ay, az = __calc_accel(p, tree, w)
    ax = post(w, p._x, vx, ax)
    ay = post(w, p._y, vx, ay)
    az = post(w, p._z, vx, az)
end

function calc_accel!(w::World{Cosmological},
    tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1},
    tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1},
    tvx::SharedArray{Float64, 1}, tvy::SharedArray{Float64, 1}, tvz::SharedArray{Float64, 1})

    w_chunks = get_chunks(w.n)
    t_chunks = get_chunks(length(tx))
    @sync for i in 1:length(workers())
        @async remotecall_wait(workers()[i], calc_accel!, w, tx,ty,tz,tax,tay,taz,tvx,tvy,tvz, w_chunks[i], t_chunks[i])
    end
end
