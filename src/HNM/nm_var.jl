using ECOS

function nm_var(m :: LCPModel,
                x :: Array{Float64,1},
                𝓘 :: Array{Int,1},  
                A :: Array{Int,1},
                E⁻ :: Array{Int,1},
                v𝓘 :: AbstractVector,
                M_𝓘𝓘 ;
                eps_active :: Float64 = 1.0e-9,
                verbose :: Bool = false)

    M = mM(m,x)
    #q = m.q
    n = size(x)
    
    d = zeros(n)
    nE⁻ = length(E⁻)
    d_E⁻ = Convex.Variable(nE⁻)
    
    AA = M_𝓘𝓘 \ Matrix(M[𝓘,E⁻])

    p = minimize(sumsquares(d_E⁻) + sumsquares( AA*d_E⁻ + v𝓘 ))

    Q = AA'*AA + I
    c = 2.0*AA'*v𝓘
    
    add_constraint!(p,  d_E⁻ >= - x[E⁻])
    AC = M[E⁻,E⁻] - M[E⁻,𝓘] * AA
    nAC = norm(AC)
    
    NAC = AC/nAC
    yx = y(m, x)  #  Éviter de recalculer
    droite = - yx[E⁻] + M[E⁻,A]*x[A] + M[E⁻,𝓘] * v𝓘

    add_constraint!(p,  NAC*d_E⁻ >= droite / nAC)

    verbose && println("ECOSSolver start")
    
    verbose && println("nE⁻ = ",nE⁻)

    solve!(p, silent=true, ECOS.Optimizer)

    verbose && println("\n\nECOSSolver stop")
    
    if any(isnan.(d_E⁻.value))
        @warn " Unfeasible QP"
        if  verbose
            println("nE⁻ = ",nE⁻)
            
            println("d_E⁻.value = ",d_E⁻.value)
            println("AA = ",AA,"\n AC = ",AC,"\n v𝓘 = ", v𝓘,
                    "\n droite = ", droite, "\n x[E⁻] = ", x[E⁻] )
            println(" sum AA² = ", sum(AA .* AA)," sum AA*v𝓘 = ",sum(AA .* v𝓘),
                    " sum v𝓘² = ", sum(v𝓘 .* v𝓘))
            println(" sum AA² = ", AA'*AA," sum AA*v𝓘 = ", AA' * v𝓘,
                    " sum v𝓘² = ", sum(v𝓘 .* v𝓘))
            
            
            unfeasibleQP = true
            @warn " Unfeasible QP"
            error(" Debug infeasible QP ")
            d[𝓘] = - v𝓘
            
            d[A] = -x[A]
            
            d[E⁻] = zeros(nE⁻)
        end
    end

    if (nE⁻ == 1)   #  annoying behaviour in Julia forces those two cases
        d[E⁻] = vec([d_E⁻.value])
    else
        d[E⁻] = vec(d_E⁻.value)
    end
    
    d[𝓘] = - (M_𝓘𝓘 \ (M[𝓘,E⁻] * d[E⁻]) ) - v𝓘

    d[A] = -x[A]
    
    return d, true, true, nE⁻
    
end
