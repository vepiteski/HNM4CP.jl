function nm_hybrid(m :: LCPModel{T},
                   x :: Vector,
                   yx :: Vector,
                   τ :: T,
                   I :: Array{Int,1},
                   A :: Array{Int,1},
                   E⁻ :: Array{Int,1},
                   I₀ :: Array{Int,1},
                   A₀ :: Array{Int,1},
                   E  :: Array{Int,1},
                   vI :: AbstractVector,
                   M_II; # :: Base.LinAlg.LU{Float64,Array{Float64,2}};
                   tol :: T = eps(),
                   τ₀ :: T = 0.01,  # Same as in linesearch!
                   eps_active :: T = 1.0e-9,
                   verbose :: Bool = false) where T
    # Newton-min-var-hybrid direction

    M = m.M

    Θ₀ = Θ(x,yx)
    
    Verbose = false
    
    n = size(x)
    
    d = T.(zeros(n))
    Θ′d = T(0.0)
    OK = true
    
    Fp = T.(zeros(n))
    d[A] = - x[A]
    
    nE⁻ = length(E⁻)
    
    # Compute the hybrid direction by trying out first to put
    # indices of E⁻ in A
    d[E⁻] = -x[E⁻]
    d[I] = - (M_II \ (M[I,E⁻]*d[E⁻]) ) - vI
    
    # try the optimality of the next point
    xp = x + d

    Md = Mtimesv(m, d)
    yp = yx + Md   #update avoids a matrix-vector product yp = y(m,xp)
    Θₓₚ = Θ(xp, yp)
    Θ′d = dir_der!(Md,d,A₀,I₀,E,x,yx,Fp)

    #verbose && println(" in Hybrid Θ(xp) = ",Θₓₚ)
    QP = false
    sizeQP = 0
    if Θₓₚ < tol # 
        x = xp
        yx = yp
        E, I₀, A₀, E⁻,E⁺ = splitsets(m, x, yx=yx, eps_active = eps_active)
        verbose && println("Solved! in hybrid")
        solved = true   #  not really use since not returned
        #  but the main algorithm will stop since θ is small
        # Test if unit stepsize satisfies Armijo condition... should receive τ₀
    elseif  ((Θₓₚ - Θ₀) <  Θ′d*τ₀) #  where     Θ₀ = Θ(x, yx)
        # again, this d will be used in the main algorithm, avoiding QP
    else
        # Compute the directional derivative estimate Θ'(x;d)
        Θ′d2 = dir_der2!(m, Md, d, τ, x, yx, eps_active = eps_active)  # formula 2.19 in the paper
        @assert !isnan(Θ′d2)
        verbose && println("Directional derivative Θ′d =",Θ′d, " Θ(m,x) = ",Θ(m,x),
                           " approx Directional derivative Θ′d2 =", Θ′d2)
        
        #if  isempty(E⁻) || ( (Θ′d2 < -0.25*Θ(m, x)) )  # JPD  d ou d2 
        if  isempty(E⁻) || ( (Θ′d2 < -0.25*Θ(x, yx)) )  # JPD  d ou d2 
            verbose && println("descent, no QP")
            Verbose && println(" d = ", d)
        else    # in last resort, solve the convex quadratic program to compute d
            d, QP, OK, sQP = nm_var(m, x, I, A, E⁻, vI, M_II, verbose=verbose)
            Md = Mtimesv(m, d)
            Θ′d = dir_der!(Md, d, A₀, I₀, E, x,yx, Fp)

            yp = yx + Md   #update avoids a matrix-vector product yp = y(m,xp)
            Θₓₚ = Θ(xp, yp)
            sizeQP = max(sizeQP, sQP)
        end

    end

    return d, Θₓₚ, Θ′d, Md, QP, OK, sizeQP
end
