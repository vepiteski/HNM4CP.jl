using Random

include("sprandsymG.jl")

export lcprandom

"""
lcpmodel(M,q,x₀), x1] = lcprandom(n,na,ni, sB, μ, L, scx, shx, scy, shy)
    
     Return data for a random linear complementarity problem 0 <= x _|_
     (M*x+q) >= 0, with a solution x and y=M*x+q satisfying
    
        n  = length(x),
        na = length(find(x < y)),
        ni = length(find(x > y)).
    
     Then, there are ne=n-na-ni indices i such that x(i)=y(i)=0. The vector
     x1 is a possible initial point for an iterative method.

    The solution is such that x[A] = scx*rand + shx and y[A] = scy*rand + shy

"""
function lcprandom(n::Int,
                   na::Int,
                   ni::Int;
                   scB :: Float64 = 500.0,
                   μ :: Float64 =  0.00005,
                   L :: Float64 = 20020.0,
                   scx::Float64 = 1.0,
                   shx::Float64 = 0.0,
                   scy::Float64 = 1.0,
                   shy::Float64 = 0.0,
                   density :: Float64 = 0.001,
                   rng :: AbstractRNG = MersenneTwister(1)
                   )
    ne = n-na-ni;

    # TODO : validate input parameters, which all need be >0 and 0<= scB <=1 
    
    #  fprintf('\nn = %i, na = %i, ne = %i, ni = %i\n',n,na,ne,ni);
    if (n <= 0) | (na < 0) | (ne < 0) | (ni < 0)
        #    fprintf('Inconsistent dimensions\n\n');
        M  = [];
        q  = [];
        x₀ = [];
        x  = [];
    end
    
    # Compute M
    
    mu = μ

    # valeurs propres aléatoires
    #v = (rand(1,n));  // diagonale contrôle cond
    #v = mu + (v-min(v))*(L-mu)/(max(v)-min(v));
    
    # valeurs propres déterministes.
    v = mu:(L-mu)/(n-1):L        
    
    # Split density in 3, Q, B and B'  (fair enough approximation)
    Q = sprandsym(n, density/3, v ;  rng = rng)

    B = sprand(rng, n,n,density/3)
    B = B - B'
    
    M = (1-scB)*Q + scB*B
    


    # generate a random solution x
    x            = zeros(n)
    x[na+ne+1:n] = scx * rand(rng, n-na-ne,1) .+ shx .+ 0.001

    Au = Vector(1:na)
    Eu = Vector(na+1:na+ne)
    Iu = Vector(na+ne+1:n)
    all = Vector(1:n)
    
    # Compute q
    
    z = M*x;

    q = zeros(n);
    
    q[1:na] = -z[1:na] + scy * rand(rng, na,1) .+ shy .+ 0.001
    q[na+1:n]       = -z[na+1:n]

    # permute randomly the components A, I and E
    rp = randperm(rng, n)
    irp = invperm(rp)
    
    x = x[rp]
    q = q[rp]
    M = M[rp,rp]
    
    # Set intial x
    
    x₀ = 10.0 * (rand(rng, n) .- 0.5)  

    As = irp[Au]
    Is = irp[Iu]
    Es = irp[Eu]
    
    return LCPModel(M, q, x₀=x₀), x, As, Is, Es
end

#  http://onlinelibrary.wiley.com/doi/10.1002/pamm.201310272/full
