function sprandsym(n,dens,eigv; rng = Xoshiro() )
    # ignore dens rounded to 1/2 or 1 for the moment
    #v = mu:(L-mu)/(n-1):L
    
    #Rseed = 590;

    szblk = Int(round(n*dens,RoundUp))
    nbblk = Int(round(n/szblk,RoundUp))

    Q = spzeros(n,n)
    for blk = 0:(nbblk-1)
        In = 1+szblk*blk : min(n,(blk+1)*szblk)
        Mi = rand(rng, length(In),length(In))
        Oi, bidon = qr(Mi) # espérons que O est aléatoire orthogonale
        Qi = Oi'*Matrix(Diagonal(eigv[In]))*Oi
        Q[In,In] = Qi
    end

    rp = randperm(rng, n)
    
    return sparse(Q[rp,rp])
end

using Test
function testsprandsym()
    n=100
    dens=0.3
    eigv=Array(1.0:100.0)

    Q = sprandsym(n,dens,eigv)

    @test isapprox(nnz(Q)/(n^2), dens, rtol = 0.1)
    @test cond(Array(Q)) ≈ 100.0

    spy(Q)
end
