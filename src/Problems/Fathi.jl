function Fathi(n :: Int = 20; T = Float64)

    MA=2.0*ones(T, n,n)
    M=triu(MA,1) + I

    M = M'*M

    q = vec(-ones(T, n))
    x = zeros(T, n)
    xs = zeros(T, n)
    xs[1] = T(1.0)

    return LCPModel(M, q, x₀ = x), xs

end
