using HNM4CP
using Test
using Printf
using SparseArrays
using LinearAlgebra

include("Problems/lcprand.jl")
include("Problems/Fathi.jl")
include("Problems/Murty.jl")
include("Problems/csizmadia.jl")
include("Problems/bg2012.jl")


@testset "HNM4CP.jl" begin
    println("\n\n Test on bg2012")
    println("------------------")
    include("nm_Test_bg2012.jl")

    println("\n\n\n Test on Csizmadia")
    println("------------------")
    include("nm_Test_Csizmadia.jl")

    println("\n\n\n Test on Murty")
    println("------------------")
    include("nm_Test_Murty.jl")

    println("\n\n\n Test on Fathi")
    println("------------------")
    include("nm_Test_Fathi.jl")

    println("\n\n\n Test on lcprand")
    println("------------------")
    include("nm_Test_lcprand.jl")
end
