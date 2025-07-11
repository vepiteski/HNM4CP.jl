using IterTools

function ispmat(M)
    nc,nr = size(M)
    if nc != nr return false
    end
    for I = subsets(1:nr)
        if (det(M[I,I]) <= 0.0)
            return false
        end
    end
    return true
end

function isp0mat(M)
    nc,nr = size(M)
    if nc != nr return false
    end
    for I = subsets(1:nr)
        if (det(M[I,I]) < 0.0)
            return false
        end
    end
    return true
end
