module HNM4CP

# Write your package code here.

using LinearAlgebra
using SparseArrays
using Complementarity
using Printf
using Suppressor
using BandedMatrices


include("HNM/types.jl")
include("HNM/nm_algo.jl")

include("Problems/lcprand.jl")
include("Problems/Fathi.jl")

export lcprandom

include("Utilities/splitSetsMod.jl")
include("Utilities/ispmat.jl")
include("Utilities/scale.jl")



end
