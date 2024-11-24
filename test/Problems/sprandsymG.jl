function sprandsym(n, dens, eigv; T = Float64, rng = Xoshiro() )
    # ignore dens rounded to 1/2 or 1 for the moment
    #v = mu:(L-mu)/(n-1):L
    
    #Rseed = 590;

    szblk = Int(round(n*dens,RoundUp))
    nbblk = Int(round(n/szblk,RoundUp))

    Q = spzeros(T, n,n)
    for blk = 0:(nbblk-1)
        In = 1+szblk*blk : min(n,(blk+1)*szblk)
        Mi = rand(rng, T, length(In),length(In))
        Oi, bidon = qr(Mi) # espérons que O est aléatoire orthogonale
        Qi = Oi'*Matrix(Diagonal(eigv[In]))*Oi
        Q[In,In] = Qi
    end

    rp = randperm(rng, n)
    
    return sparse(Q[rp,rp])
end

using Test
function testsprandsym(;T=Float64)
    n=50
    dens=0.3
    eigv=Array{T}(1.0:50.0)

    Q = sprandsym(n,dens,eigv, T = T)

    @test isapprox(nnz(Q)/(n^2), dens, rtol = 0.1)
    @test cond(Array(Q)) ≈ 50.0

    spy(Q)
end
