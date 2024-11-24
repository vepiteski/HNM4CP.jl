export bg2012, BG2012

function bg2012(n,verb;T)

##
# [M,q,x,x1] = bg2012(n,verb)
#
# Return data for the BenGharbia-Gilbert problem with n variables (n
# muts be ≥ 3). On return, M is a sparse P-matrix. Data are taken from
#
#   I. Ben Gharbia, J.Ch. Gilbert (2012), "Nonconvergence of the plain
#   Newton-min algorithm for linear complementarity problems with a
#   P-matrix, Mathematical Programming, 134, 349-364,
#   http://dx.doi.org/10.1007/s10107-010-0439-6.
#
# If verb is false, no printing.

  M  = T.([])	# matrix M of the problem
  q  = T.([])	# vector q of the problem
  x  = T.([])	# solution to the problem
  x1 = T.([])	# suggested initial point

  if isempty(n) | (n < 3)
      if verb;
          @printf("\n\n### bg2012: the dimension n = %i must be ≥ 3\n\n",n);
      end
    return
  end

# Introductory message

  #cl = clock;
  #dstr = sprintf('%i-%i-%i, %i:%i:%i',cl(3),cl(2),cl(1),cl(4),cl(5),fix(cl(6)));
  if verb
    @printf("\nbg2012 problem ");
    @printf("\n. n = %i\n",n);
  end

# Compute M

  M = spzeros(n,n) + I

    if isodd(n)   	# n is odd
        alpha  = T(2)
        a      = alpha*ones(n-1)
        M += spdiagm(-1 => a)
        M[1,n] = alpha
    else			# n is even
        alpha  = T(4)/T(3)
        a      = alpha*ones(n-1)
        M += spdiagm(-1 => a)
        M[1,n] = alpha

        beta  = T(big"0.5")
        b  = beta*ones(n-2)
        M += spdiagm(-2 => b)
        b  = beta*ones(2)
        M += spdiagm(n-2 => b)
    end
    
# Set q

  q = ones(n);

# Set x (the solution)

  x = zeros(n);

# Set x1 (initial point)

  x1    = zeros(n);
  x1[1] = T(-1);

    return M,q,x,x1
end


function BG2012(n :: Int =64; T = Float64)
    M, q, x, x1 = bg2012(n, false, T=T)

    return LCPModel(M, q, x₀ = x1), x
end


    
