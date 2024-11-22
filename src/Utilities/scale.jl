function scale(M0 :: AbstractMatrix,
               q0 :: Array{Float64,1},
               xs :: Array{Float64,1},  # solution of unscaled problelm
               x₀ :: Array{Float64,1};  # starting point of unscaled problem
               S :: Symbol = S3,
               tol :: Float64 = 0.5,
               ntype :: Real = 2)

    M = M0
    q = q0
    x = x₀
    xsS = xs
    ΣL = ones(size(xs))
    ΣC = ones(size(xs))

    # S0 is no scaling, just return the input data...
    if S == :S1
        #nM0 = norm(Array(M0), ntype)
        nM0 = norm((M0), ntype)
        M = M0 / nM0
        q = q0 / nM0
    elseif S == :S2  # scaling the lines of M and q, no effect on xs or x₀
        M, ΣL = scaleLt(M0, ΣL, ntype = ntype)
        q = q0 .* ΣL
    elseif S == :S3 # scaling both the lines and columns of M, q and x
        # very simple version where the columns are scaled, then the lines, and that's all.
        # Columns
        done = false
        M = M0
        
        while ~done
            M, ΣCi = scaleC(M, ΣC, ntype = ntype)
            ΣC = ΣC .* ΣCi
            M, ΣLi = scaleLt(M, ΣL, ntype = ntype)
            ΣL = ΣL .* ΣLi
            scaling = maximum(abs.(1.0 .- [ΣCi; ΣLi]))
            #println("         | 1 - scaling| = ",scaling)
            done = scaling < tol
        end
        x₀ = x₀ ./ ΣC
        xs = xs ./ ΣC
        q = q0 .* ΣL
    end

    return M, q, xs, x₀, ΣL, ΣC

end


function scaleC(M0 :: AbstractMatrix, ΣC :: AbstractVector; ntype :: Real = 2)
    n, = size(M0)
    Σ = ones(n)   #  diagonal scaling matrix stored as a vector
    for i=1:n
        Σ[i] /= norm(M0[:,i])
    end
    M = M0 * spdiagm(0 => Σ)
    ΣC = ΣC .* Σ

    return M, Σ
end

function scaleL(M0 :: AbstractMatrix, ΣL :: AbstractVector; ntype :: Real = 2)

    n, = size(M0)
    Σ = ones(n)   #  diagonal scaling matrix stored as a vector
    for i=1:n
        Σ[i] /= norm(M0[i,:])
    end
    M = spdiagm(0 => Σ) * M0
    ΣL = ΣL .* Σ

    return M, Σ
end

function scaleLt(M0 :: AbstractMatrix, ΣL :: AbstractVector; ntype :: Real = 2)

    n, = size(M0)
    Σ = ones(n)   #  diagonal scaling matrix stored as a vector

    Mt = (M0')
    for i=1:n
        Σ[i]  /= norm(Mt[:,i])
    end
    M = spdiagm(0 => Σ) * M0
    ΣL = ΣL .* Σ

    return M, Σ
end
