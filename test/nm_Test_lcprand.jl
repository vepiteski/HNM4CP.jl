
using HNM4CP

using Random
using SparseArrays
using LinearAlgebra
using Printf
#
pbtype = :lcprand
n = 128

na = Int(round(0.3*n))
ni = Int(round(0.5*n))
scaleB = 0.1
density = 0.1

scx = 01.0
shx = 00000.0

scy = 0001.0
shy = 0000.0

T = Float64

println("Testing lcprandom with n = ",n, ", density = ", density, ", scaleB = ", scaleB)
println(" and scale x = ",scx, ", shift x = ",shx, ", scale y = ",scy, ", shift y = ",shy )
m0, xs0, As, Is, Es = lcprandom(n, na, ni, scB = scaleB, μ = 0.001, L = 1850.0, density = density,
                                scx=scx, shx=shx, scy=scy, shy=shy)

println(" Actual density = ", nnz(m0.M)/(n^2))
println(" cond(M) = ", cond(Array(m0.M)))

E, I₀, A₀, E⁻,E⁺ = splitsets(m0, xs0)
@printf("%10s  %17s        %16.2e  %3i  %3i  %3i    \n", "True solution","",Θ(m0, xs0),length(As),length(Is),length(Es))

@printf("%10s  %19s        %9.2e  %3i  %3i  %3i  %3i   \n", "Splitsets solution","  ",Θ(m0, xs0),length(A₀),length(I₀),length(E⁺),length(E⁻))


@test issetequal(E, Es)
@test isapprox(Θ(m0,xs0),T(0.0),atol=sqrt(eps(T)))
#
m=m0


label = "HPNM"
@printf("Variant    niter   time  #M*vector Θ(m,xsol) |A₀| |I₀| |E⁺| |E⁻|    #E  #QP   MoyQP   #tval  max|E⁻|  \n")

exectime = time()
xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m)
exectime = time() - exectime

@printf("%5s     %5i   %5.2f  %8i %9.2e  %3i  %3i  %3i  %3i    %3i  %3i  %6.1f  %5i  %3i\n",
        label,niter, exectime, m.nf, Θ(m,xsol),length(A₀),length(I₀),length(E⁺),length(E⁻), nbE, nbQP, mQP, nbtval, maxE )

