# function exec(p::Particle, hist::History)
#     vx = 0.0
#     vy = 0.0
#     vz = 0.0
#     ax = 0.0
#     ay = 0.0
#     az = 0.0
#     tree = hist.tree[1]
#     ax, ay, az = calc_accel(p,tree,hist.w)
#     for i in 1:(length(hist.dt)-1)
#         # kick half dt
#         dt = hist.dt[i]
#         vx += ax*dt/2
#         vy += ay*dt/2
#         vz += az*dt/2
#         # drift dt
#         p = addxyz(p, vx*dt, vy*dt, vz*dt)
#         # kick another half dt
#         tree = hist.tree[i+1]
#         ax, ay, az = calc_accel(p,tree,hist.w)
#         vx += ax*dt/2
#         vy += ay*dt/2
#         vz += az*dt/2
#     end
#     p
# end
