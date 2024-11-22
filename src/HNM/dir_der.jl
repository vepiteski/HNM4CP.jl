function dir_der!(Md, d, AA₀, II₀, EE, x, yx, Fp)
    E, I₀, A₀ = splitsets(x, yx)
    Fp[A₀] = d[A₀]
    Fp[E] = min.(d[E],Md[E])
    Fp[I₀] = Md[I₀]
    return F(x,yx) ⋅ Fp
end



function dir_der0!(Md, d, x, yx)
    E, I₀, A₀ = splitsets(x, yx)
    Fp = zeros(size(x))
    Fp[A₀] = d[A₀]
    Fp[E] = min.(d[E],Md[E])
    Fp[I₀] = Md[I₀]
    return F(x,yx) ⋅ Fp
end

# Formula 2.19b and ρ given by formula 2.20
function dir_der2!(m:: LCPModel{T}, Md, d, τ, x, yx; eps_active = T(1e-9) ) where T
    
    E, II, A, E⁻,E⁺ = splitsets(m, x, τ, yx = yx, eps_active = eps_active )

    I₀ = II[findall(abs.(yx[II]) .< eps_active)]
    A₀ = A[findall(abs.(x[A]) .< eps_active)]

    ρxd = zeros(size(x))
    I⁺ = setdiff(II, I₀)
    A⁺ = setdiff((A ∪ E⁺), A₀)


    ρxd[A⁺] = (x[A⁺] + d[A⁺]) ./ x[A⁺]
    ρxd[I⁺] = (yx[I⁺] + Md[I⁺]) ./ yx[I⁺]
    ρxd[E⁻] = max.((x[E⁻] + d[E⁻]) ./ x[E⁻], (yx[E⁻] + Md[E⁻]) ./ yx[E⁻])

    replace!(ρxd, NaN => T(0.0))
    replace!(ρxd, Inf => T(0.0))
    replace!(ρxd, -Inf => T(0.0))

    D =  (T(1.0) .- ρxd)
    # estimate of Θ'  formula 2.19a in the paper
    Θ′d2 = -sum(D .* F(x,yx).^2)

    return Θ′d2
end
