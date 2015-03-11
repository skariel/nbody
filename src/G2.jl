# simple Gadget2 save/load functionality

function save(w::nb.World, fn::String; t=0.0)
    open(fn, "w") do f

        #
        #      HEADER!
        #
        ########################

        write(f, [UInt(w.n), 0, 0, 0, 0, 0])
        write(f, [w.particles[1]._m, 0.0, 0.0, 0.0, 0.0, 0.0])
        write(f, t)
        write(f, 1./t-1.)
        write(f, 0) # flag SFR
        write(f, 0) # flag feedback
        write(f, [UInt(w.n), 0, 0, 0, 0, 0]) # Nall
        write(f, 0) # flag cooling
        write(f, 1) # number of files in each snapshot
        write(f, 0.0) # box size
        write(f, 1.0) # Ω_0
        write(f, 0.0) # Ω_Λ
        write(f, 2./3) # h
        # padding with zeros
        write(f, [UInt8(0) for i in position(f):255])


        #
        #      POSITIONS!
        #
        ########################
        for p in w.particles
            write(f, [Float32(p._x), Float32(p._y), Float32(p._z)])
        end

        #
        #      VELOCITIES!
        #
        ########################
        for i in 1:w.n
            write(f, [Float32(w.vx[i]), Float32(w.vy[i]), Float32(w.vz[i])])
        end

        #
        #      IDs!
        #
        ########################
        for i in 1:w.n
            write(f, i)
        end

    end
end

function load(fn::String)
    open(fn, "r") do f
        n = [Int(read(f,UInt)) for i in 1:6]
        m = [read(f,Float64) for i in 1:6]
        t = read(f,Float64)
        z = read(f,Float64)
        fsfr = read(f,Int)
        ffbk = read(f,Int)
        n = [Int(read(f,UInt)) for i in 1:6]
        fcool = read(f,Int) # flag cooling
        nfiles = read(f,Int) # number of files in each snapshot
        bsize = read(f,Float64) # box size
        Ω0 = read(f,Float64)
        ΩΛ = read(f,Float64)
        h = read(f,Float64)
        # padding...
        [read(f, UInt8) for i in position(f):255 ]

        # positions...
        pos_x = Float64[]
        pos_y = Float64[]
        pos_z = Float64[]
        for i in 1:w.n
            push!(pos_x, Float64(read(f,Float32)))
            push!(pos_y, Float64(read(f,Float32)))
            push!(pos_z, Float64(read(f,Float32)))
        end

        # velocities
        vel_x = Float64[]
        vel_y = Float64[]
        vel_z = Float64[]
        for i in 1:w.n
            push!(vel_x, Float64(read(f,Float32)))
            push!(vel_y, Float64(read(f,Float32)))
            push!(vel_z, Float64(read(f,Float32)))
        end

        id = [Int(read(f,Int)) for i in 1:w.n]
        @show id[11]
    end
end
