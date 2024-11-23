export Murty

function Murty(n :: Int = 20;
               variant :: Symbol = :Z, T = Float64)

    MA=2.0*ones(T,n,n)
    M=tril(MA,-1) + I

    x=zeros(T, n)
    q = vec(-ones(T, n,1))

    xs=zeros(T, n)
    xs[1] = 1.0

    if variant != :Z
        # According to Murty original analysis,
        #https://public.websites.umich.edu/~murty/books/linear_complementarity_webbook/
        # other variants
        # but since q is huge, useful only for n≤66, otherwise
        # PATH crashes badly and exits Julia!
        q = vec(zeros(T, n,1))
        q[1] = -2.0^n
        for i=2:n
            q[i] = q[i-1] - 2.0^(n-(i-1))
        end
        if variant == :p1
            x = vec(ones(T, n,1))
        else # assume that variant == :n1
            x = -vec(ones(T, n,1))
        end
        xs[1] = T(2.0^n)
    end

    return LCPModel(M, q, x₀ = x), xs
    
end
