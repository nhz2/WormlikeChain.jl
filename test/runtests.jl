using WormlikeChain
using Test
using RNGTest
using HypothesisTests
using Distributions
using Suppressor



"""
Roughly tests if the 1d array x is normally distributed iid. 
    Not super rigorous
"""
function testifnormal(x)
    sets= [x[:],
            x[1:2:end],
            x[2:2:end],
            √0.5(x[1:2:end] .+ x[2:2:end]),
            x[1:4:end],
            x[2:4:end],
            x[3:4:end],
            x[4:4:end],
            √0.25(x[1:4:end] .+ x[2:4:end] .+ x[3:4:end] .+ x[4:4:end]) ]
    for a in sets
        t= 0.0
        @suppress_err begin
            t=ApproximateOneSampleKSTest(a,Normal())
        end
        @test pvalue(t) > 0.005
    end
end

@testset "WormlikeChain.jl" begin
    # Write your tests here.
end

@testset "random_utils.jl" begin
    N= 1_000_000
    x= zeros(2*N)
    key= UInt64(123215)
    ctr1= UInt64(1)
    ctr2= UInt64(0) 
    #generate a bunch of normal random numbers
    for i in 1:N
        x[2i], x[2i - 1] = WormlikeChain.randn_2x64(key, ctr1, ctr2+i)
    end
    testifnormal(x)
    for i in 1:N
        x[2i], x[2i - 1] = WormlikeChain.randn_2x64(key, ctr1+i, ctr2)
    end
    testifnormal(x)

    x= zeros(Float32,4*N)
    for i in 1:N
        x[4i], x[4i-1], x[4i-2], x[4i-3] = WormlikeChain.randn_4x32(key, ctr1, ctr2+i)
    end
    testifnormal(x)
    for i in 1:N
        x[4i], x[4i-1], x[4i-2], x[4i-3] = WormlikeChain.randn_4x32(key, ctr1+i, ctr2)
    end
    testifnormal(x)
    
end