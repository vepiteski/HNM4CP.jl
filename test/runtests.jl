using HNM4CP
using Test

@testset "HNM4CP.jl" begin
    println("Test on bg2012")
    println("------------------")
    include("nm_Test_bg2012.jl")

    println("Test on Csizmadia")
    println("------------------")
    include("nm_Test_Csizmadia.jl")

    println("Test on Murty")
    println("------------------")
    include("nm_Test_Murty.jl")

    println("\n\n\n Test on fathi")
    println("------------------")
    include("nm_Test_Fathi.jl")

    println("\n\n\n Test on lcprand")
    println("------------------")
    include("nm_Test_lcprand.jl")
end
