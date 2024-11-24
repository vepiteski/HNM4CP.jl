using HNM4CP

using Printf
##########   Fathi   ##############
#

T = Float64

n=256

println("Testing Fathi with n = ",n)
m, xs = Fathi(n; T=T)

@test isapprox(Θ(m,xs),T(0.0),atol=sqrt(eps(T)))
#

@printf("Variant    niter   time  #M*vector Θ(m,xsol) |A₀| |I₀| |E⁺| |E⁻|    #E  #QP   MoyQP   #tval  max|E⁻|  \n")

label = "HPNM"

exectime = time()

xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m);

exectime = time() - exectime

@test isapprox(Θ(m,xsol),T(0.0),atol=sqrt(eps(T)))
@printf("%5s     %5i   %5.2f  %8i %9.2e  %3i  %3i  %3i  %3i    %3i  %3i  %6.1f  %5i  %3i\n",
        label,niter, exectime, m.nf, Θ(m,xsol),length(A₀),length(I₀),length(E⁺),length(E⁻), nbE, nbQP, mQP, nbtval, maxE )

@test isapprox(Θ(m,xsol),T(0.0),atol=sqrt(eps(T)))

