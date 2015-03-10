# not much work for newtonian forces other than common stuff...
# mainly forwarding

@inline calculate_accel_on_particle!(w::World{Newtonian}, i::Int64) =
    __calculate_accel_on_particle!(w,i)
@inline calc_accel!(w::World{Newtonian}, rng::UnitRange{Int64}) =
    __calc_accel!(w, rng)
@inline calc_accel!(w::World{Newtonian}, tx::SharedArray{Float64, 1}, ty::SharedArray{Float64, 1}, tz::SharedArray{Float64, 1}, tax::SharedArray{Float64, 1}, tay::SharedArray{Float64, 1}, taz::SharedArray{Float64, 1}, w_rng::UnitRange{Int64}, t_rng::UnitRange{Int64}) =
    __calc_accel!(w,tx,ty,tz,tax,tay,taz,w_rng,t_rng)
@inline calc_accel(p::Particle, tree::CompiledOctTree{Particle}, w::World{Newtonian}) =
    __calc_accel(p, tree, w)
