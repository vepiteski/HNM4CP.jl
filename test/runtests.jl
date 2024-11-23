using HNM4CP
using Test

@testset "HNM4CP.jl" begin
    # Write your tests here.
    println("Test on lcprand")
    println("------------------")
    include("nm_Test_lcprand.jl")
    println("Test on fathi")
    println("------------------")
    include(nm_Test_Fathi.jl")
end
