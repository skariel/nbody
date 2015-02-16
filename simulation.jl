function calcdt(w::World)
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

function kick(w::World, dt::Float64)
    @inbounds for i in 1:w.n
        w.vx[i] += w.ax[i]*dt
        w.vy[i] += w.ay[i]*dt
        w.vz[i] += w.az[i]*dt
    end
end

function drift(w::World, dt::Float64)
    @inbounds for i in 1:w.n
        w.x[i] += w.vx[i]*dt
        w.y[i] += w.vy[i]*dt
        w.z[i] += w.vz[i]*dt
    end
end
