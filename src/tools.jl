"""
$(SIGNATURES)

Calculate unit normal vector

2D
"""
function unit_normal(p1::T, p2::T) where {T<:AbstractVector}
    Δ = p2 .- p1
    l = norm(Δ) + 1e-6

    return [-Δ[2], Δ[1]] ./ l
end

"""
$(SIGNATURES)

3D
"""
function unit_normal(p1::T, p2::T, p3::T) where {T<:AbstractVector}
    v1 = p2 .- p1
    v2 = p3 .- p1

    n = cross(v1, v2)
    l = norm(n) + 1e-6

    return n ./ l
end
