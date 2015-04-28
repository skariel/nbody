Ha(a, w) = H0*sqrt(w.Ω0/a/a/a+w.ΩΛ)

FD(a1,a2,w) = quadgk(a->1./Ha(a,w)/a/a/a,a1,a2)[1]
FK(a1,a2,w) = quadgk(a->1./Ha(a,w)/a/a,a1,a2)[1]

# for convenience...

FAC1(w) = 0.5*H0*H0*w.Ω0
