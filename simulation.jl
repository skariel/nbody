function calc_dt(w::World)
    mindt2 = 1.e30 # infinity, ha!
    @inbounds for i in 1:w.n
        const a2 = w.ax[i]*w.ax[i] + w.ay[i]*w.ay[i]
        const dyn_dt2 = sqrt(w.smth2/a2)*w.dtfrac*w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    sqrt(mindt2)
end

function kick!(w::World; dt=0.0)
    @inbounds for i in 1:w.n
        w.vx[i] += w.ax[i]*dt
        w.vy[i] += w.ay[i]*dt
        w.vz[i] += w.az[i]*dt
    end
end

function drift!(w::World; dt=0.0)
    @inbounds for i in 1:w.n
        const dx = w.vx[i]*dt
        const dy = w.vy[i]*dt
        const dz = w.vz[i]*dt
        w.particles[i] = withxyz(w.particles[i], dx, dy, dz)
    end
end

function exec!(sim::Simulation)
    reset!(sim)
    calc_accel!(sim.w)
    break_time = false
    tic()
    while true
        sim.step += 1
        sim.dt = calc_dt(sim.w)
        if !sim.limit_by_steps
            if sim.t+sim.dt > sim.tf
                sim.dt = sim.tf - sim.t
                break_time = true
            end
        end
        print("s=",sim.step," t=",sim.t," dt=",sim.dt)
        tic()

        kick!(sim.w, dt=sim.dt/2)
        drift!(sim.w, dt=sim.dt)
        kick!(sim.w, dt=sim.dt/2)
        calc_accel!(sim.w)

        sim.t += sim.dt
        elapsed = toq()
        print(" /",elapsed,"s\n")

        sim.limit_by_steps && sim.step >= sim.stepf && break
        break_time && break
    end
    println("\n--- Done!\n")
    toc()
    nothing
end
