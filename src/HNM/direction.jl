function direction(m :: LCPModel{T},
                   x :: Array{T,1},             # current iterate
                   yx :: Array{T,1};            # M*x + q
                   variant :: Symbol = :hybrid, # HNM or PNM
                   τ :: T = 1.0e-9,             # tolerance for indices in E⁻
                   τ₀:: T = 0.01,               # Armijo parameter,          
                                                #  used in the hybrid direction to bypass QP
                   eps_active :: T = 1.0e-9,
                   tol :: T = 1e-9,             # stopping criterion,
                                                #  used in hybrid to bypass further computations
                   verbose ::Int = 0,
                   EinI :: Bool = false) where T


    #@assert isapprox(yx, y(m,x))
    # Index sets at x
    E, I₀, A₀, E⁻,E⁺ = splitsets(m, x, T(τ), yx = yx, eps_active = eps_active)

    #
    # Choice: put all E⁺ in A for now
    Î = I₀

    # useless if-then-else but compensates for inneficient set operations
    # in Julia

    if isempty(E⁺) A = A₀   else   A = A₀ ∪ E⁺   end
    #

    #
    # try to put E in I
    if EinI
        #@show "E in I"
        Î = I₀ ∪ E⁺
        A = A₀
    end
    #

    M = mM(m,x)
    
    # Factorization
    #
    # trouble with banded matrices since M[I,I] seems to
    # convert to a temporary full matrix :(
    if isa(M,BandedMatrix)
        println("banded")
        l = bandwidth(M,1)
        u = bandwidth(M,2)
        sI = length(Î)
        MII = BandedMatrix(BandedMatrices.Zeros(sI,sI),(l,u))
        MII.data[:,:] = M.data[:,Î]
        M_II = MII
        #M_II = factorize(MII)
    else
    #
        M_II = M[Î,Î]
        M_II = factorize(M_II)
    end

    d  = zeros(T,size(x))
    Md = zeros(T,size(x))
    Θxd = T(0.0)
    Θ′d = T(0.0)
    QP = false
    sQP = 0
    OK = true
    
    Fp = zeros(T,size(x))
    d[A] = - x[A]

    start = time()
    vI = M_II \ (yx[Î] - M[Î,A]*x[A])
    nE⁻ = length(E⁻)
    if isempty(E⁻)
        d[Î] = -vI
        Md = Mtimesv!(m, d, Md)
        Θ′d = dir_der!(Md, d, A₀, I₀, E, x, yx, Fp)
        yx +=  Md
        x += d
        Θxd = Θ(x,yx)
    else 
        if variant == :var
            d, QP, OK, sQP = nm_var(m, x, Î, A, E⁻, vI, M_II, verbose=verbose>2)
            Md = Mtimesv(m, d)
            Θ′d = dir_der!(Md, d, A₀, I₀, E, x, yx, Fp)
            if (Θ′d > 0.0) && (verbose != 0)
                @printf("%9.2e  %3i  %3i  %3i  %3i   \n",
                        Θ′d, length(A₀),length(I₀),length(E⁺),length(E⁻))
                
                (verbose==0) || error( " :var Ascent direction")
            end

            Θxd = Θ(x+d, yx + Md)  # save a md  Θxd = Θ(m, x+d)
        elseif variant == :hybrid
            d, Θxd, Θ′d, Md, QP, OK, sQP = nm_hybrid(m, x, yx, τ, Î, A, E⁻, I₀, A₀, E, vI, M_II,
                                                     verbose=verbose>1, τ₀ = τ₀, tol = tol, eps_active = eps_active)
        else    # for now, only two variants, :var and :hybrid
            error("unimplemented direction")
        end
    end
    if (Θ′d > 0.0) && (verbose != 0)
        @printf("%9.2e  %3i  %3i  %3i  %3i   \n",
                Θxd,length(A₀),length(I₀),length(E⁺),length(E⁻))
        
        @warn( "Ascent direction") ; @show  Θ′d 
    end
    
    return d, Θxd, Θ′d, Md, QP, OK, A₀, I₀, E⁺, E⁻, E, sQP
end
