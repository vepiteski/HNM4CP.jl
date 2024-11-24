# HNM4CP

<p align="center">
<img src="./work-progress_5578703.png" width=200 height=200>
</p>

<a href="https://www.freepik.com/icons/work-in-progress">Icon by Freepik</a>

[![Build Status](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package will provide a basket of tools to deal with Complementarity Problems. Now, it is limited to Linear Complementarity Problems in the so called standard form

$$0\le x \perp Mx+q \ge 0$$

**Abstract data type**


The complementarity problems are encapsulated in a data type, all complying to ``` AbstractCPModel`{T}``. For the moment, the only concrete type is for standard LCP, ```LCPModel{T} <: AbstractCPModel{T}``` available in the file ```src/HNM/types.jl```. The T refers to floating point precision which may be selected as Float32, Float64, Float128, BigFloat, Double64 or even other future complying float representations.

**PNM and HNM solvers**

The solvers implement the algorithms in DFG2024 https://hal.science/hal-02306526 . PNM is a Newton-min algorithm made globally convergent using Polyhedral projections. The hybrid HNM version uses tests to avoid as much as possible the polyhedral computations for the sake of efficiency, while maintaining the global convergence property.

**Examples of complementarity problems**

In src/Problems we provide a collection of random generator and problems from applications which are used to test the solvers.

**Interfaces to other solvers**

We also provide interfaces to compare HNM to other available solvers.

------------------------------

**Installation and usage**

For now, the package is in construction and is not registered in the Julia ecosystem. Therefore, it is somewhat more complicated to install and use. 


In the REPL, type  "]"  to enter the Package manager subshell.
Then, type in

      add https://github.com/vepiteski/HNM4CP.jl

or if you intend to modify and contribute,

      dev https://github.com/vepiteski/HNM4CP.jl

The package was built using Julia 1.11.1 but hopefully will work using other (recent) versions.

Still in the Pkg subshell, type in

      test HNM4CP

this will trigger tests, hopefully successful.

**Example usage**

Inspect the folder Test to imitate the test scripts to use the HNM solver.

For example, copy any test, say `nm_Test_Fathi.jl` in a fresh folder. Start Julia in this folder.

    * In the julia REPL, type in `using Test`.

    * Then, again in the Julia REPL, typein `include("nm_Test_Fathi.jl")`, which should perform 
    as expected.

You are now in position to modify the nm_Test_Fathi.jl script.


Much more and much clearer instructions are in preparation!!!


