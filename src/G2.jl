# simple Gadget2 save/load functionality

function save(w::World, fn::String; t=0.0)
    open(fn, "w") do f

        #
        #      HEADER!
        #
        ########################

        write(f, Int32(256))
        write(f, [UInt32(0), UInt32(w.n), UInt32(0), UInt32(0), UInt32(0), UInt32(0)])
        write(f, [0.0, w.particles[1]._m, 0.0, 0.0, 0.0, 0.0])
        write(f, t)
        write(f, 1./t-1.)
        write(f, Int32(0)) # flag SFR
        write(f, Int32(0)) # flag feedback
        write(f, [Int32(0), Int32(w.n), Int32(0), Int32(0), Int32(0), Int32(0)]) # Nall
        write(f, Int32(0)) # flag cooling
        write(f, Int32(1)) # number of files in each snapshot
        write(f, 0.0) # box size
        write(f, 1.0) # Ω_0
        write(f, 0.0) # Ω_Λ
        write(f, 2./3) # h
        # padding with zeros
        while position(f) < 256+4
            write(f, UInt8(0x00))
        end
        @show position(f)
        write(f, Int32(256))

        #
        #      POSITIONS!
        #
        ########################
        write(f, Int32(4*w.n*3))
        for p in w.particles
            write(f, [Float32(p._x), Float32(p._y), Float32(p._z)])
        end
        write(f, Int32(4*w.n*3))

        #
        #      VELOCITIES!
        #
        ########################
        write(f, Int32(4*w.n*3))
        for i in 1:w.n
            write(f, [Float32(w.vx[i]), Float32(w.vy[i]), Float32(w.vz[i])])
        end
        write(f, Int32(4*w.n*3))

        #
        #      IDs!
        #
        ########################
        write(f, Int32(4*w.n))
        for i in 1:w.n
            write(f, Int32(i))
        end
        write(f, Int32(4*w.n))

    end
end

function load(fn::String)
    open(fn, "r") do f
        ff=read(f,Int32)
        n = [Int(read(f,UInt32)) for i in 1:6]
        m = [read(f,Float64) for i in 1:6]
        t = read(f,Float64)
        z = read(f,Float64)
        fsfr = read(f,Int32)
        ffbk = read(f,Int32)
        n = [Int(read(f,Int32)) for i in 1:6]
        fcool = read(f,Int32) # flag cooling
        nfiles = read(f,Int32) # number of files in each snapshot
        bsize = read(f,Float64) # box size
        Ω0 = read(f,Float64)
        @show ΩΛ = read(f,Float64)
        @show h = read(f,Float64)
        # padding...
        while position(f) < 256+4
            read(f, UInt8)
        end
        @assert ff==read(f,Int32)

        # positions...
        ff=read(f,Int32)
        pos_x = Float64[]
        pos_y = Float64[]
        pos_z = Float64[]
        for i in 1:sum(n)
            push!(pos_x, Float64(read(f,Float32)))
            push!(pos_y, Float64(read(f,Float32)))
            push!(pos_z, Float64(read(f,Float32)))
        end
        @assert ff==read(f,Int32)

        # velocities
        ff=read(f,Int32)
        vel_x = Float64[]
        vel_y = Float64[]
        vel_z = Float64[]
        for i in 1:sum(n)
            push!(vel_x, Float64(read(f,Float32)))
            push!(vel_y, Float64(read(f,Float32)))
            push!(vel_z, Float64(read(f,Float32)))
        end
        @assert ff==read(f,Int32)

        ff=read(f,Int32)
        id = [Int(read(f,Int32)) for i in 1:sum(n)]
        @assert ff==read(f,Int32)
        @show id[13]
        # TODO: build world here...
        smth=0.2
        opening_alpha=0.5
        dtfrac=0.35
        space=Cosmological
        return World([Particle(pos_x[i],pos_y[i],pos_z[i],maximum(m)) for i in 1:sum(n)], smth, opening_alpha, dtfrac, space)
    end
end
