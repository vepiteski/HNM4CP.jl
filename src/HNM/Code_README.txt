Solvers for LCPs.

The nm_var implements PNM and HNM solvers. Its main component are organized
as follows.

The main algorithm is in nm_algo.jl. The algorithm calls direction.jl, itself
forking in PNM (nm_var.jl) and HNM (nm_hybrid.jl). Then linesearch.jl is called.

dir_der.jl include computations of directional derivative, an exact one and
a bound corresponding to equation (2.19) in the paper.

JPD October, 14 2024.
