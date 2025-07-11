# HNM4CP

<p align="center">
<img src="./work-progress_5578703.png" width=200 height=200>
</p>

<a href="https://www.freepik.com/icons/work-in-progress">Icon by Freepik</a>

[![Build Status](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package will provide a basket of tools to deal with Complementarity Problems. Now, it is limited to Linear Complementarity Problems in the so called standard form,

$$0\le x \perp (Mx+q) \ge 0.$$



**Abstract data type**


The complementarity problems are encapsulated in a data type, all complying to ``` AbstractCPModel{T}```. For the moment, the only concrete type is for standard LCP, ```LCPModel{T} <: AbstractCPModel{T}``` available in the file ```src/HNM/types.jl```. The T refers to floating point precision which may be selected as Float32, Float64, Float128, BigFloat, Double64 or even other future complying float representations.

**PNM and HNM solvers**

The solvers implement the algorithms in DFG2024 https://hal.science/hal-02306526 . PNM is a Newton-min algorithm made globally convergent using Polyhedral projections. The hybrid HNM version uses tests to avoid as much as possible the polyhedral computations for the sake of efficiency, while maintaining the global convergence property. The main solver is in ```/src/HNM/nm_algo.jl```. There is a ```/src/NMM/Code_README.txt``` briefly describing the solver.

**Examples of complementarity problems**

In ```test/Problems``` we provide a collection of random generator and problems from applications which are used to test the solvers.

**Interfaces to other solvers**

We will also provide soon interfaces to compare HNM to other available solvers.

------------------------------

**Installation and usage**

For now, the package is in construction and is not registered in the Julia ecosystem. Therefore, it is somewhat more complicated to install and use. 


In the REPL, type  "]"  to enter the Package manager subshell.
Then, type in

      add https://github.com/vepiteski/HNM4CP.jl

The package was built using Julia 1.11.1 but hopefully will work using other (recent) versions.

Still in the Pkg subshell, type in

      test HNM4CP

this will trigger tests, hopefully successful.

**Example usage**

**A simple example**

Let a LCP be defined by
```
M = [ 0.11  22.5 ;
     -0.85   2.5 ]
q = [-45.6;
     -2   ]
```
The solver receives a problem produced using the function 
```
m = LCPModel(M,q)
```
or, if specifying an initial point is desired,
```
x = [-0.2; 
      2.019]
m = LCPModel(M,q,x₀=x)
```

The simplest application consists in calling 
```
julia> nm_algo(m)
([3.556701030927835, 2.009278350515464], Int64[], [1, 2], Int64[], Int64[], 3, 0, 0, 0, 0, 0, 0)
```
A more useful call will allow to interpret the result:
```
julia> xsol, A₀, I₀, E⁺, E⁻, niter, nbE, nbQP, nbtval, maxE, sQP, mQP  = nm_algo(m);

julia> xsol
2-element Vector{Float64}:
 3.556701030927835
 2.009278350515464
```
which yields the solution point. If only the solution point is desired, a shorter call works.
```
(xsol, )= nm_algo(m);
```


If an execution trace is desired, the parameter verbose may be specified.
- verbose = 0    default
- verbose = 1    prints options summary
- verbose = 2    prints iterations
- verbose = 3    prints also line search iterations

For example, verbose = 2 gives

```
julia> nm_algo(m, verbose = 2);
Newton min variants algorithm.

 number of variables:             2
 maximum iterations:              4000
 direction variant:               hybrid
 Put indices of E⁺ in I :         false
 kink help to line search:        first
 non monotone line search memory: 10
 eps_active:                      1.0e-9 Float64
 absolute stopping tolerance:     2.220446049250313e-16 Float64

Niter   Θ         |A₀| |I₀|  |E⁺| |E⁻|   stepsize   QP       dimQP
    0   2.06e+00    2    0    0    0  
    1   2.06e+00    2    0    0    0    1.2207e-04  false     0
    2   1.09e-01    1    1    0    0    7.8125e-03  false     0
    3   0.00e+00    0    2    0    0    1.0000e+00  false     0
Solved! 
    3   0.00e+00    0    2    0    0   
```





**Example from the tests**

Inspect the folder ```test``` to imitate the test scripts to use the HNM solver.

For example, copy any test, say `nm_Test_Fathi.jl` in a fresh folder. Start Julia in this folder.

    * In the julia REPL, type in `using Test`.

    * Then, again in the Julia REPL, typein `include("nm_Test_Fathi.jl")`, which should perform 
    as expected.

You are now in position to modify the nm_Test_Fathi.jl script.


Much more and much clearer instructions are in preparation!!!


