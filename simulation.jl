function calc_dt(w::World)
    mindt2 = 1.e30 # infinity, ha!
    @inline for i in 1:w.n
        const a2 = w.ax[i]*w.ax[i] + w.ay[i]*w.ay[i]
        const dyn_dt2 = sqrt(w.smth2/a2)*w.dtfrac*w.dtfrac
        if dyn_dt2 < mindt2
            mindt2 = dyn_dt2
        end
    end
    sqrt(mindt2)
end

function kick!(w::World; dt=0.0)
    @inline for i in 1:w.n
        w.vx[i] += w.ax[i]*dt
        w.vy[i] += w.ay[i]*dt
        w.vz[i] += w.az[i]*dt
    end
end

function drift!(w::World; dt=0.0)
    @inline for i in 1:w.n
        const dx = w.vx[i]*dt
        const dy = w.vy[i]*dt
        const dz = w.vz[i]*dt
        w.particles[i] = addxyz(w.particles[i], dx, dy, dz)
    end
end

# TODO: make this parallel!
function calc_accel!(sim::Simulation)
    buildtree!(sim.w, sim.tree)
    calc_accel!(sim.w)
end

function exec!(sim::Simulation; use_brute_force=false, silent=false, accumulate_history=false)
    reset!(sim)
    tic()
    calc_accel!(sim)
    break_time = false
    hist = History(CompiledOctTree{Particle}[], Float64[], sim.w)
    while true
        sim.step += 1
        sim.dt = calc_dt(sim.w)
        if accumulate_history
            push!(hist.tree, deepcopy(sim.w.tree))
            push!(hist.dt, sim.dt)
        end
        if !sim.limit_by_steps
            if sim.t+sim.dt > sim.tf
                sim.dt = sim.tf - sim.t
                break_time = true
            end
        end
        !silent && print("s=",sim.step," t=",sim.t," dt=",sim.dt)
        tic()

        kick!(sim.w, dt=sim.dt/2)
        drift!(sim.w, dt=sim.dt)
        if !use_brute_force
            calc_accel!(sim)
        else
            calc_accel_brute_force!(sim)
        end
        kick!(sim.w, dt=sim.dt/2)

        sim.t += sim.dt
        elapsed = toq()
        !silent && print(" /",elapsed,"s\n")

        sim.limit_by_steps && sim.step >= sim.stepf && break
        break_time && break
    end
    if accumulate_history
        push!(hist.tree, deepcopy(sim.w.tree))
        push!(hist.dt, sim.dt)
    end
    !silent && println("\n--- Done!\n")
    silent? toq() : toc()
    hist
end
