# simple Gadget2 save/load functionality

function save_G2(w::World, fn::String, t)
    open(fn, "w") do f

        #
        #      HEADER!
        #
        ########################

        write(f, int32(256))
        write(f, [uint32(0), uint32(w.n), uint32(0), uint32(0), uint32(0), uint32(0)])
        write(f, [0.0, w.particles[1]._m, 0.0, 0.0, 0.0, 0.0])
        write(f, t)
        write(f, 1./t-1.)
        write(f, int32(0)) # flag SFR
        write(f, int32(0)) # flag feedback
        write(f, [int32(0), int32(w.n), int32(0), int32(0), int32(0), int32(0)]) # Nall
        write(f, int32(0)) # flag cooling
        write(f, int32(1)) # number of files in each snapshot
        write(f, 0.0) # box size
        write(f, w.Ω0) # Ω_0
        write(f, w.ΩΛ) # Ω_Λ
        write(f, 2./3) # h
        # padding with zeros
        while position(f) < 256+4
            write(f, uint8(0x00))
        end
        write(f, int32(256))

        #
        #      POSITIONS!
        #
        ########################
        write(f, int32(4*w.n*3))
        for p in w.particles
            write(f, [float32(p._x), float32(p._y), float32(p._z)])
        end
        write(f, int32(4*w.n*3))

        #
        #      VELOCITIES!
        #
        ########################
        write(f, int32(4*w.n*3))
        for i in 1:w.n
            write(f, [
                float32(w.vx[i]/(t.^(3/2))),
                float32(w.vy[i]/(t.^(3/2))),
                float32(w.vz[i]/(t.^(3/2)))
            ])
        end
        write(f, int32(4*w.n*3))

        #
        #      IDs!
        #
        ########################
        write(f, int32(4*w.n))
        for i in 1:w.n
            write(f, int32(i))
        end
        write(f, int32(4*w.n))

    end
end

function load_G2(fn::String)
    open(fn, "r") do f
        ff=read(f,Int32)
        n = [int(read(f,Uint32)) for i in 1:6]
        m = [read(f,Float64) for i in 1:6]
        t = read(f,Float64)
        z = read(f,Float64)
        fsfr = read(f,Int32)
        ffbk = read(f,Int32)
        n = [int(read(f,Int32)) for i in 1:6]
        fcool = read(f,Int32) # flag cooling
        nfiles = read(f,Int32) # number of files in each snapshot
        bsize = read(f,Float64) # box size
        Ω0 = read(f,Float64)
        ΩΛ = read(f,Float64)
        h = read(f,Float64)
        # padding...
        while position(f) < 256+4
            read(f, Uint8)
        end
        @assert ff==read(f,Int32)

        # positions...
        ff=read(f,Int32)
        pos_x = Float64[]
        pos_y = Float64[]
        pos_z = Float64[]
        for i in 1:sum(n)
            push!(pos_x, float64(read(f,Float32)))
            push!(pos_y, float64(read(f,Float32)))
            push!(pos_z, float64(read(f,Float32)))
        end
        @assert ff==read(f,Int32)

        # velocities
        ff=read(f,Int32)
        vel_x = Float64[]
        vel_y = Float64[]
        vel_z = Float64[]
        for i in 1:sum(n)
            push!(vel_x, float64(read(f,Float32)))
            push!(vel_y, float64(read(f,Float32)))
            push!(vel_z, float64(read(f,Float32)))
        end
        @assert ff==read(f,Int32)

        ff=read(f,Int32)
        id = [int(read(f,Int32)) for i in 1:sum(n)]
        @assert ff==read(f,Int32)
        # TODO: build world here...
        smth=0.2
        opening_alpha=0.5
        dtfrac=0.35
        space=Cosmological
        w = World([Particle(pos_x[i],pos_y[i],pos_z[i],maximum(m)) for i in 1:sum(n)], smth, opening_alpha, dtfrac, space, Ω0, ΩΛ)
        for i in 1:w.n
            w.vx[i] = vel_x[i].*(t.^(3/2))
            w.vy[i] = vel_y[i].*(t.^(3/2))
            w.vz[i] = vel_z[i].*(t.^(3/2))
        end
        return w
    end
end
