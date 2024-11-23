using HNM4CP
using Test

@testset "HNM4CP.jl" begin
    # Write your tests here.
    println("Test on Murty")
    println("------------------")
    include("nm_Test_Murty.jl")

    println("Test on fathi")
    println("------------------")
    include("nm_Test_Fathi.jl")

    println("Test on lcprand")
    println("------------------")
    include("nm_Test_lcprand.jl")
end
