
function estimDiff(M)
    n, = size(M)
    l₀ = Inf
    for i=1:n^2
        z = (rand(n) .- 0.5)
        z = z / norm(z,Inf)
        zMz = z .* (M*z)
        l₀ = min(l₀, sqrt( maximum(zMz) ) )
    end
    return l₀
end

function estimDiff2(M,p)
    n, = size(M)
    βₚ = 0.0
    for i=1:n^2
        d = rand(n)
        D = Diagonal(d)
        A = (I - D + D*M)\D
        
        βₚ = max(βₚ, norm(A,p) )
    end
    return βₚ*norm(M,p)
end

function estimDiff3(M,p)
    n, = size(M)
    Kₚ = 0.0
    for i=1:n^2
        d = rand(n)
        D = Diagonal(d)
        A = (I - D + D*M)
        
        Kₚ = max(Kₚ, norm(A,p)*norm(inv(Matrix(A)),p) )
    end
    return Kₚ
end

