function linesearch!(m :: LCPModel{T},
                     x :: Vector,
                     Θₓ :: T,
                     d :: Vector,
                     Md :: Vector,
                     yx :: Vector,
                     Θxd :: T,
                     Θ′d :: T,
                     valΘ :: Vector,
                     niter :: Int;
                     kinkStep :: Symbol = :first,
                     N :: Int = 1,
                     verbose :: Bool = false,
                     τ₀ :: T = 0.01) where T

    #################################
    #  line search
    #################################
    # simple Armijo backtracking
    #               (non monotone)
    # Θ′d is assumed valid for d
    #################################
    verbose && println(" Line Search ")
    t = T(1.0)
    Θₓ = Θ(x, yx)
    valΘ[mod(niter,N) + 1] = Θₓ
    Θ₀ = maximum(valΘ)
    xp = x + t*d
    nbktr = 0
    nbtval = 0
    yp = y(m, xp)
    Θp = Θ(xp, yp)

    # éviter les calculs impliquant m dans les prints
    verbose && println("     Θ′d = ",Θ′d, " Θ(x) = ",Θ(x,yx), " Θ(xp) = ", Θp, " Θxd = ", Θxd)
    verbose && println("     Θ₀ = ", Θ₀, "    Θₓ = ", Θₓ)
    bktrmax = 70; fact = T(big"0.5")   #  for 0.5, 50 is  enough

    Θₓ = Θ(x, yx)

    # Uncomment for systematic unit stepsize
    #return xp, Θp, t, 1, nbktr, true
    verbose && @printf("         Θ(x⁺)        t        Θ′d*t*τ₀ \n")
    verbose && @printf("      %9.2e   %9.2e   %9.2e\n", Θp, t, Θ′d*t*τ₀)
    
    while  ((Θp - Θ₀) >  Θ′d*t*τ₀) && (nbktr<bktrmax)
        t *= fact
        xp = x + t*d
        yp = yx + t*Mdnone
        Θp = Θ(xp, yp)
        
        #verbose && println("         Θ(xp) = ", Θp, "  t = ",t, " Θ0 = ",Θ₀, "Θ′d*t*τ₀ = ",Θ′d*t*τ₀)
        
        verbose && @printf("      %9.2e   %9.2e   %9.2e\n", Θp, t, Θ′d*t*τ₀)
        nbktr += 1
    end

    if nbktr == bktrmax
        verbose && println(" Max iter: Θ(xp) = ", Θp, "  t = ",t, " Θ0 = ",Θ₀, "Θ′d*t*τ₀ = ",Θ′d*t*τ₀)
        
        return xp, Θp, t, 0, nbktr, false
    end

    # Kinks computation
    kink = (x-yx) ./ (Md-d) #  éviter de recalculer y(x)  et  M*d
    kink_pos = kink[findall(kink .> T(0.0))]
    kink_pos = kink_pos[findall(isfinite.(kink_pos))]
    if !isempty(kink_pos)
        t_min = minimum(kink_pos)
    else
        t_min = T(0.0)
    end

   

    # test to go at least as far as the first kink
    if (kinkStep == :first)
        if (t < t_min) & ( Θ₀  >=  Θ(x+t_min*d, yx+t_min*Md) - Θ′d*t_min*τ₀ )
            t = t_min
            xp = x + t*d
            nbtval += 1
        end
    end
    # test to go at least as far as the best kink
    if kinkStep == :best 
        if !isempty(kink_pos)
            val_kink = zeros(size(kink_pos))
            for i=1:length(kink_pos)
                xi = x + kink_pos[i] * d
                yi = yx + kink_pos[i] * Md
                val_kink[i] = Θ(xi, yi)
            end
            val_min, ind_min = findmin(val_kink)
            tval_min = kink_pos[ind_min]
        else
            tval_min = 0.0
        end
        if (t < tval_min) && ( Θ(m, x+t*d) > val_min )
            t = tval_min
            xp = x + t*d
            nbtval += 1
        end
    end

    yp = yx + t*Md
    θp = Θ(xp, yp)
    verbose && println(" At end: Θ(x⁺) = ", Θp, "  t = ", t )

    return xp, Θp, t, nbtval, nbktr, true
end
