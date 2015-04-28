# filtering stuff

function filter_r{T<:SpaceType}(w::World{T}, r::Float64)
    particle_arr = Particle[]
    vx = Float64[]
    vy = Float64[]
    vz = Float64[]
    rt_2 = r*r
    xs = [p._x for p in w.particles]
    d = 0.5*(maximum(xs) + minimum(xs))
    for i in 1:w.n
        p = w.particles[i]
        x = p._x - d
        y = p._y - d
        z = p._z - d
        r2 = x*x+y*y+z*z
        if r2 < rt_2
            push!(particle_arr, Particle(x, y, z, p._m))
            push!(vx, w.vx[i])
            push!(vy, w.vy[i])
            push!(vz, w.vz[i])
        end
    end
    mvx = mean(vx)
    mvy = mean(vy)
    mvz = mean(vz)
    @show smth = sqrt(w.smth2)
    @show opening_alpha = sqrt(w.opening_alpha2)
    @show dtfrac = w.dtfrac
    @show space = typeof(w.space)
    @show Ω0 = w.Ω0
    @show ΩΛ = w.ΩΛ
    new_world = World(particle_arr, smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)
    for i in 1:new_world.n
        new_world.vx[i] = vx[i]-mvx
        new_world.vy[i] = vy[i]-mvy
        new_world.vz[i] = vz[i]-mvz
    end
    new_world
end

function filter_rate{T<:SpaceType}(w::World{T}, rate::Int64)
    particle_arr = Particle[]
    vx = Float64[]
    vy = Float64[]
    vz = Float64[]
    for i in 1:rate:w.n
        p = w.particles[i]
        push!(particle_arr, p)
        push!(vx, w.vx[i])
        push!(vy, w.vy[i])
        push!(vz, w.vz[i])
    end
    # fixing mass
    for i in 1:length(particle_arr)
        p = particle_arr[i]
        particle_arr[i] = Particle(p._x,p._y,p._z,p._m*w.n/length(particle_arr))
    end
    mvx = mean(vx)
    mvy = mean(vy)
    mvz = mean(vz)
    @show smth = sqrt(w.smth2)
    @show opening_alpha = sqrt(w.opening_alpha2)
    @show dtfrac = w.dtfrac
    @show space = typeof(w.space)
    @show Ω0 = w.Ω0
    @show ΩΛ = w.ΩΛ
    new_world = World(particle_arr, smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)
    for i in 1:new_world.n
        new_world.vx[i] = vx[i]-mvx
        new_world.vy[i] = vy[i]-mvy
        new_world.vz[i] = vz[i]-mvz
    end
    new_world
end

# redshift space

function get_z_coords(x::Float64, y::Float64, z::Float64, vx::Float64, vy::Float64, vz::Float64)
  r2 = x*x+y*y+z*z+0.001
  r = sqrt(r2)
  rx = x/r
  ry = y/r
  rz = z/r
  vr = vx*rx+vy*ry+vz*rz
  zx = x + vr*rx / H0
  zy = y + vr*ry / H0
  zz = z + vr*rz / H0
  zx, zy, zz
end

function get_c_coords(x::Float64, y::Float64, z::Float64, vx::Float64, vy::Float64, vz::Float64)
  r2 = x*x+y*y+z*z+0.001
  r = sqrt(r2)
  rx = x/r
  ry = y/r
  rz = z/r
  vr = vx*rx+vy*ry+vz*rz
  zx = x - vr*rx / H0
  zy = y - vr*ry / H0
  zz = z - vr*rz / H0
  zx, zy, zz
end

function to_z_space!(w::World)
    for i in 1:w.n
        p = w.particles[i]
        zx,zy,zz = get_z_coords(p._x,p._y,p._z,w.vx[i],w.vy[i],w.vz[i])
        w.particles[i] = Particle(zx,zy,zz,p._m)
    end
end

function to_c_space!(w::World)
    for i in 1:w.n
        p = w.particles[i]
        zx,zy,zz = get_c_coords(p._x,p._y,p._z,w.vx[i],w.vy[i],w.vz[i])
        w.particles[i] = Particle(zx,zy,zz,p._m)
    end
end
