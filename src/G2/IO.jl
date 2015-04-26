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

        N = w.n
        tmp_all = zeros(Float32, 3N)
        for i in 1:N
            tmp_all[(i-1)*3+1] = float32(w.particles[i]._x)
            tmp_all[(i-1)*3+2] = float32(w.particles[i]._y)
            tmp_all[(i-1)*3+3] = float32(w.particles[i]._z)
        end
        write(f, tmp_all)
        write(f, int32(4*w.n*3))

        #
        #      VELOCITIES!
        #
        ########################
        write(f, int32(4*w.n*3))
        for i in 1:w.n
            tmp_all[(i-1)*3+1] = float32(w.vx[i]/(t.^(3/2)))
            tmp_all[(i-1)*3+2] = float32(w.vy[i]/(t.^(3/2)))
            tmp_all[(i-1)*3+3] = float32(w.vz[i]/(t.^(3/2)))
        end
        write(f, tmp_all)
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
    nothing
end

function load_G2(fn::String)
    open(fn, "r") do f
        ff=read(f,Int32)
        n = [int(read(f,Uint32)) for i in 1:6]
        @show n
        m = [read(f,Float64) for i in 1:6]
        @show m
        @show t = read(f,Float64)
        @show z = read(f,Float64)
        @show z
        @show t
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

        N = sum(n)
        tmp_all = zeros(Float32, 3N)
        read!(f, tmp_all)

        for i in 1:sum(n)
            push!(pos_x, tmp_all[(i-1)*3+1])
            push!(pos_y, tmp_all[(i-1)*3+2])
            push!(pos_z, tmp_all[(i-1)*3+3])
        end
        @assert ff==read(f,Int32)

        # velocities
        ff=read(f,Int32)
        vel_x = Float64[]
        vel_y = Float64[]
        vel_z = Float64[]
        read!(f, tmp_all)
        for i in 1:sum(n)
            push!(vel_x, tmp_all[(i-1)*3+1])
            push!(vel_y, tmp_all[(i-1)*3+2])
            push!(vel_z, tmp_all[(i-1)*3+3])
        end
        @assert ff==read(f,Int32)

        #ff=read(f,Int32)
        #id = [int(read(f,Int32)) for i in 1:sum(n)]
        #@assert ff==read(f,Int32)

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
