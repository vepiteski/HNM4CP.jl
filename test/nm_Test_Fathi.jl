using HNM4CP

using Printf
##########   Fathi   ##############
#

#for n=128:256

T = Float64

n=256 #Itérations: n jusqu'à 166, 114 de 167 à 
println("Testing Fathi with n = ",n)
m, xs = Fathi(n; T=T)    
@test isapprox(Θ(m,xs),T(0.0),atol=sqrt(eps(T)))
#
verbose = 0
tol=T(BigFloat("1e-30"))
@printf("Variant           niter   time  #M*vector Θ(m,xsol) |A₀| |I₀| |E⁺| |E⁻|     #E  #QP   MoyQP   #tval  max|E⁻|  \n")

label = "HPNM"
#for eps_test in [1.0e-5,1.0e-6,1.0e-7,1.0e-8,1.0e-9,1.0e-10,1.0e-11,1.0e-12]
eps_test = 1.0e-7

exectime = time()
xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m, verbose = verbose, kink = :first, τ₀ = T(0.0001), eps_active = T(eps_test), τ = T(1.0e-7), tol=tol, maxiter = 2*n, EinI = true);
#xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m, verbose = 0, kink = :first, τ₀ = 0.0001, eps_active = eps_test, τ = 1.0e-7);
exectime = time() - exectime
@test isapprox(Θ(m,xsol),T(0.0),atol=sqrt(eps(T)))

@printf("%5s %6.1e     %5i   %5.2f  %8i %9.2e  %3i  %3i  %3i  %3i    %3i  %3i  %6.1f  %5i  %3i\n",
        label, eps_test,niter, exectime, m.nf, Θ(m,xsol),length(A₀),length(I₀),length(E⁺),length(E⁻), nbE, nbQP, mQP, nbtval, maxE )
@test isapprox(Θ(m,xsol),T(0.0),atol=sqrt(eps(T)))
#end

