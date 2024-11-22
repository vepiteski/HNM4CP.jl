export splitsets

"""
    splitsets in details...
"""
function splitsets(x :: Vector, τ :: T, yx :: Vector ; eps_active :: T = 1e-9) where T
    n = length(x)
    xmyx = x - yx
    Eτ  = findall( abs.(xmyx) .< max(τ, eps_active) )
    
    # delicate treatement of near zero xmyx[i]
    E₀  = findall(abs.(xmyx) .< eps_active)
    I   = findall(xmyx  .> eps_active)   # no use of Identity here, simpler to keep the notation I
    A   = findall(xmyx  .< -eps_active)
    E   = setdiff([1:n], A ∪ I)
    # I, A and E form a partition of [1:n]
        
    I⁻  = I[findall(x[I]  .< 0)]
    I⁺  = setdiff(I,I⁻)

    A⁻  = A[findall(yx[A] .< 0)]
    A⁺  = setdiff(A,A⁻)

    AI⁻ = I⁻ ∪ A⁻
    Eτ⁻ = Eτ[findall(max.(x[Eτ],yx[Eτ]) .< 0.0)]
    AI⁺ = I⁺ ∪ A⁺
    E₀⁺ = E₀[findall(max.(x[E₀],yx[E₀]) .>= 0.0)]

    S = Eτ
    UI = setdiff(I,Eτ⁻)
    UA = setdiff(A,Eτ⁻)
    
    return Eτ, UI, UA, Eτ⁻, E₀⁺, I, A, I⁻, A⁻
end


#
function splitsets(m:: LCPModel{T}, x :: Vector, τ :: T ; yx :: Vector = y(m,x), eps_active :: T = 1e-9)    where T        

    return splitsets(x, τ, yx, eps_active = eps_active)
end



function splitsets(m:: LCPModel, x :: Vector; yx :: Vector = y(m,x), eps_active :: T = 1.0e-9) where T

    return splitsets(x, yx, eps_active = eps_active)
end    
#

function splitsets(x :: Vector, yx :: Vector; eps_active :: T = 1.0e-9) where T
    # TODO:  rewrite by calling the above with τ=0
    n = length(x)
    xmyx = x - yx
    E  = findall(abs.(xmyx) .<= eps_active)
    I₀ = findall(xmyx .>  eps_active);
    A₀ = findall(xmyx .< -eps_active);

    E⁻ = E[findall(x[E] .< 0.0)]
    E⁺ = setdiff(E, E⁻)

    return E, I₀, A₀, E⁻,E⁺
end

#=
function test()
    M = rand(5,5)
    q = rand(5)

    x = rand(5) .- 0.5
    
    y = M*x + q
    y[1] = x[1] + 1e-10
    
    E1, I1, A1, E1⁻,E1⁺ = splitsets(x, y)
    E2, I2, A2, E2⁻,E2⁺, I, A, I⁻, A⁻= splitsets(x, 0.0, y)
    println("test1",E1, E2⁺)
    @test E1 == E2
    println("test2",A1, A)
    @test A1 == A2
    println("test3",I1, I)
    @test I1 == I2
end
=#
