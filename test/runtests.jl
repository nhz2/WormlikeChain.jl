using WormlikeChain
using Test
using RNGTest
using HypothesisTests
using Distributions

@testset "WormlikeChain.jl" begin
    # Write your tests here.
end

@testset "random_utils.jl" begin
    N= 10_000_000
    x= zeros(2*N)
    key= UInt64(123213)
    ctr1= UInt64(1)
    ctr2= UInt64(0) 
    #generate a bunch of normal random numbers
    for i in 1:N
        x[2i], x[2i - 1] = WormlikeChain.randn_2x64(key, ctr1, ctr2+i)
    end
    t=ApproximateOneSampleKSTest(x[2:2:end],Normal())
    @test pvalue(t) > 0.2
    t=ApproximateOneSampleKSTest(x[1:2:end],Normal())
    @test pvalue(t) > 0.2
    t=ApproximateOneSampleKSTest(√0.5(x[1:2:end] .+ x[2:2:end]),Normal())
    @test pvalue(t) > 0.2

end