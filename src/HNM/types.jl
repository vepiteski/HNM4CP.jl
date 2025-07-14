"""
Ancestor of all complementarity models
"""
abstract type AbstractCPModel{T} end;

"""
Linear Complementarity Problem

Place holder for LCP of the form 0 ≤ x ⟂ (Mx + q) ≥ 0

Three access functions (also inplace! results):

        y(m,x)       = M*x + q
        Mtimesv(m,v) = M*v
        mM(m)        = M

        reset!(m)    resets the counters
"""
mutable struct LCPModel{T} <: AbstractCPModel{T}
    M :: AbstractMatrix{T}
    q :: AbstractVector{T}

    x₀ :: AbstractVector{T}
    nf :: Int  # Counts the number of M*v products, the only O(n²) operation
    ng :: Int  # Not useful for LCPs since M is constant
end

function LCPModel( M :: AbstractMatrix{T},
                   q :: AbstractVector{T};
                   x₀ :: AbstractVector{T} = zeros(T,size(q)) ) where T
    return LCPModel(M, q, x₀, 0, 0)
end

" resets the counters "
function reset!(m::LCPModel)
    m.nf = 0
    m.ng = 0
end


# LCP functions
function y(m::LCPModel, x::AbstractVector)
    m.nf += 1
    return m.M*x + m.q
end

function y!(m::LCPModel, x::AbstractVector, yx::AbstractVector)
    m.nf += 1
    yx = m.M*x + m.q
    return yx
end

function Mtimesv(m::LCPModel, v::AbstractVector)
    m.nf += 1
    res :: AbstractVector = m.M*v
    return res
end

function Mtimesv!(m::LCPModel, v::AbstractVector, res::AbstractVector)
    m.nf += 1
    res =  m.M*v
    return res
end



function mM(m::LCPModel, x::AbstractVector)
    m.ng += 1
    return m.M
end

F(x::AbstractVector,y::AbstractVector)   = min.(x,y)
F(m::LCPModel,      x::AbstractVector)   = F(x,y(m,x))

function Θ(x::AbstractVector{T},y::AbstractVector{T}) where T
    return T(0.5) * F(x,y)⋅F(x,y)
end

Θ(m::LCPModel,      x::AbstractVector)   = Θ(x,y(m,x))

export AbstractCPModel, LCPModel, F, Θ, y, y!, nM, Mtimesv!, Mtimesv, reset!
