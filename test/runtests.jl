using WormlikeChain
using Test
using HypothesisTests
using Distributions
using Suppressor
using LinearAlgebra



"""
Roughly tests if the 1d array x is normally distributed iid. 
    Not super rigorous, but hopefully catches some issues
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

include("test_interface.jl")

include("test_cpukernels.jl")

@testset "bonded_forces_utils.jl" begin
    @test 1 == WormlikeChain.next_beadid(1,(1,2))
    @test 2 == WormlikeChain.next_beadid(1,(1,10))
    @test 4 == WormlikeChain.next_beadid(3,(1,10))
    @test 1 == WormlikeChain.next_beadid(9,(1,10))
    @test 1 == WormlikeChain.next_beadid(9,(1,10,15))
    @test 11 == WormlikeChain.next_beadid(10,(1,10,15))
    @test 10 == WormlikeChain.next_beadid(14,(1,10,15))

    @test 1 == WormlikeChain.prev_beadid(1,(1,2))
    @test 9 == WormlikeChain.prev_beadid(1,(1,10))
    @test 2 == WormlikeChain.prev_beadid(3,(1,10))
    @test 8 == WormlikeChain.prev_beadid(9,(1,10))
    @test 9 == WormlikeChain.prev_beadid(1,(1,10,15))
    @test 14 == WormlikeChain.prev_beadid(10,(1,10,15))
    @test 13 == WormlikeChain.prev_beadid(14,(1,10,15))

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
        x[4i-3:4i] .= WormlikeChain.randn_4x32(key, ctr1, ctr2+i)
    end
    testifnormal(x)
    for i in 1:N
        x[4i-3:4i] .= WormlikeChain.randn_4x32(key, ctr1+i, ctr2)
    end
    testifnormal(x)

    x= zeros(Float32,6*N)
    for i in 1:N
        #make sure to inc ctr2 by at least 2 because each 4 normal samples uses up 1 ctr2
        x[6i-5:6i] .= WormlikeChain.randn_32(key, ctr1, ctr2+i*2,Val(6))
    end
    testifnormal(x)
    for i in 1:N
        x[6i-5:6i] .= WormlikeChain.randn_32(key, ctr1+i, ctr2,Val(6))
    end
    testifnormal(x)
    
end