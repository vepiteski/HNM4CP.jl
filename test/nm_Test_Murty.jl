using HNM4CP

using Printf

##########   Murty   ##############
#
n=64

T=Float64
tol=T(BigFloat("1e-30"))

println("Testing Murty with n = ",n)
m, xs = Murty(n, variant = :Z, T=T)
@test isapprox(Θ(m,xs),T(0.0),atol=sqrt(eps(T)))


eps_test = 1e-7
verbose = 1
#

label = "HPNM"

@printf("Variant           niter   time  #M*vector Θ(m,xsol) |A₀| |I₀| |E⁺| |E⁻|     #E  #QP   MoyQP   #tval  max|E⁻|  \n")

exectime = time()
xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m, verbose = verbose, kink = :first, τ₀ = T(0.0001), eps_active = T(eps_test), τ = T(1.0e-7), tol=tol, maxiter = Int(6*n), EinI = false)
exectime = time() - exectime
@test isapprox(Θ(m,xsol),T(0.0),atol=sqrt(eps(T)))

@printf("%10s      %5i   %5.2f  %8i %9.2e  %3i  %3i  %3i  %3i    %3i  %3i  %6.1f  %5i  %3i\n",
        label, niter, exectime, m.nf, Θ(m,xsol),length(A₀),length(I₀),length(E⁺),length(E⁻), nbE, nbQP, mQP, nbtval, maxE )


