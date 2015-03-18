type World{T<:SpaceType}
    tree::CompiledOctTree{Particle}
    particles::SharedArray{Particle, 1}
    vx::SharedArray{Float64, 1}
    vy::SharedArray{Float64, 1}
    vz::SharedArray{Float64, 1}
    ax::SharedArray{Float64, 1}
    ay::SharedArray{Float64, 1}
    az::SharedArray{Float64, 1}
    dt::SharedArray{Float64, 1}
    n::Int64
    opening_alpha2::Float64
    smth2::Float64
    dtfrac::Float64
    space::T
    mindt::Float64
end
