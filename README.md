# HNM4CP
![image info](./work-in-progress.png =100x20)

[![Build Status](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vepiteski/HNM4CP.jl/actions/workflows/CI.yml?query=branch%3Amain)

In construction. Will provide a Hybrid Newton-Min solver for complementarity problems. Next step, provide the solver for LCPs.

First test incorporated.

For now, no easy way to use this package.

In the REPL, type  "]"  to enter the Package manager subshell.
Then, type in

      dev https://github.com/vepiteski/HNM4CP.jl

This should install the package. Still in the Pkg subshell, type in

      test HNM4CP

this should trigger tests, hopefully successful.

Then, inspect the folder Test to imitate the test scripts to use the HNM solver.

For example, copy any test, say `nm_Test_Fathi.jl` in a fresh folder. Start Julia in this folder.

    * In the julia REPL, type in `using Test`.

    * Then, again in the Julia REPL, typein `include("nm_Test_Fathi.jl")`, which should perform as expected.

You are now in position to modify the nm_Test_Fathi.jl script.




Much more and much clearer instructions are in preparation!!!


