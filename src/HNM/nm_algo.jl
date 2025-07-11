include("direction.jl")
include("nm_var.jl")
include("nm_hybrid.jl")
include("linesearch.jl")
include("dir_der.jl")
using Convex
#using SCS
using Base.Iterators

export nm_algo


function nm_algo(m :: LCPModel{T};
                 tol :: AbstractFloat = eps(T),
                 maxiter :: Int = 4000,
                 dir :: Symbol = :hybrid,    # HNM or PNM
                 kink :: Symbol = :first,    # Kink to aim in the LS
                 N :: Int = 10,               # Memory for non monotone LS
                 verbose :: Int = 0,
                 τ = T(1.0e-7),      # Tolerance for E⁻
                 #τ :: AbstractFloat = T(1.0e-7),      # Tolerance for E⁻
                 τ₀:: AbstractFloat = T(0.0001),        # Armijo descent parameter
                 eps_active  = T(1.0e-9),
                 EinI :: Bool = false) where T
    # Newton min general algorithm

    M = mM(m,m.x₀)
    q = m.q
    x = m.x₀

    (verbose >0) && println("Newton min variants algorithm.\n")
    (verbose >0) && println(" number of variables:             ", length(x))
    (verbose >0) && println(" maximum iterations:              ", maxiter)
    (verbose >0) && println(" direction variant:               ", dir)
    (verbose >0) && println(" Put indices of E⁺ in I :         ", EinI)
    (verbose >0) && println(" kink help to line search:        ", kink)
    (verbose >0) && println(" non monotone line search memory: ", N)
    (verbose >0) && println(" eps_active:                      ", eps_active, " ",typeof(eps_active))
    (verbose >0) && println(" absolute stopping tolerance:     ", tol, " ",typeof(tol),"\n\n")
    
    Θ′d = 0.0   # place holder for the directional derivative of Θ
    
    niter = 0
    solved = false
    unfeasibleQP = false
    

    #   a few counters
    nbQP = 0    # number of QP required
    sQP  = 0    # size of QP
    kQP = 0
    nbE = 0     # number of iterations with non     empty E
    nbtval = 0  # number of iterations where the Armijo step is replaced by
                # the first or best kink
    maxE = 0    # the maximum size |E| encountered

    valΘ = T(-Inf) * ones(T,N)
    yx = y(m, x)
    Θₓ = Θ(x, yx)

    (verbose >1) && @printf("Niter   Θ         |A₀| |I₀|  |E⁺| |E⁻|   stepsize   QP       dimQP\n")  
    E, I₀, A₀, E⁻,E⁺ = splitsets(m, x, T(τ), yx = yx, eps_active = eps_active)
    (verbose >1) && @printf("%5i  %9.2e  %3i  %3i  %3i  %3i  \n",
                                niter,Θₓ,length(A₀),length(I₀),length(E⁺),length(E⁻))
    # Julia nonsense
    Θ⁺ = Θₓ
    yx⁺ = yx
    
    Epred, Ipred, Apred = splitsets(m, x, yx = yx)

    while ~solved && (niter < maxiter)        

        d, Θxd, Θ′d, Md, QP, OK, A₀, I₀, E⁺, E⁻, E, sizeQP = direction(m, x, yx,
                                                                       variant = dir,
                                                                       τ = τ,
                                                                       eps_active = eps_active,
                                                                       tol = tol,
                                                                       verbose = verbose,
                                                                       τ₀ = τ₀,
                                                                       EinI = EinI)

        OK || ( error(" Problem in computing the direction."))
        if QP  nbQP += 1 end
        if ~isempty(E)
            nbE += 1
            maxE = max(maxE,length(E⁻))
        end
        sQP = max(sQP,sizeQP)
        kQP += sizeQP
        t = 1.0

        OK = true
        if Θxd <= tol
            x⁺, Θ⁺, yx⁺ = x+d, Θxd, yx + Md
        else
            x⁺, Θ⁺, t, tval, nbktr, OK = linesearch!(m, x, Θₓ, d, Md, yx, T(Θxd), T(Θ′d),
                                                     valΘ, niter, kinkStep = kink,
                                                     N = N, τ₀ = τ₀, verbose = verbose > 2)

            #(verbose >2) && @printf("Niter   Θ         |A₀| |I₀|  |E⁺| |E⁻|   stepsize   QP       dimQP\n")  
            #(verbose >1) && @printf("%5i  %9.2e  %3i  %3i  %3i  %3i   %11.4e  %s  %4i\n",
            #                        niter,Θₓ,length(A₀),length(I₀),length(E⁺),length(E⁻),t,
            #                        QP ? "true" : "false", sizeQP)
            
            nbtval += tval
            yx⁺ = yx + t*Md
            yx⁺ = y(m,x⁺)
        end

        niter += 1
        x, Θₓ, yx = x⁺, Θ⁺, yx⁺

        (verbose >2) && @printf("Niter   Θ         |A₀| |I₀|  |E⁺| |E⁻|   stepsize   QP       dimQP\n")  
        (verbose >1) && @printf("%5i  %9.2e  %3i  %3i  %3i  %3i   %11.4e  %s  %4i\n",
                                niter,Θₓ,length(A₀),length(I₀),length(E⁺),length(E⁻),t,
                                QP ? "true" : "false", sizeQP)
        
       
        solved = (Θₓ <= tol)
                
    end
    if niter < maxiter
        # Index sets at x
        E, I₀, A₀, E⁻,E⁺ = splitsets(m, x, yx = yx, eps_active = eps_active)

        (verbose >1) && @printf("Solved! \n")
        (verbose >1) && @printf("%5i  %9.2e  %3i  %3i  %3i  %3i   \n",
                                niter,Θₓ,length(A₀),length(I₀),length(E⁺),length(E⁻))

    else
        verbose>0 && @printf("Failure, max iteration reached\n")
    end

    if kQP==0 moyQP =0 else moyQP= kQP/nbQP end
    
    return x, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, moyQP
end
