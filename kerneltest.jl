using KernelAbstractions
using CUDA
using CUDAKernels
using Random123
using WormlikeChain

const N= 256

@kernel function foo(A)
    I = @index(Global)
    uI= I%UInt64
    @inbounds A[I]= 0f0
    sharedA = @localmem Float32 N
    @inbounds sharedA[I] = 0
    for i in 1:1000000
        ns= WormlikeChain.randn_4x32(UInt64(0), uI, i%UInt64)
        @inbounds sharedA[I]+= ns[1]
    end
    @inbounds A[I]= sharedA[I]
end

fookh = foo(CPU(), 8)
fookd = foo(CUDADevice(), N)

A_h = ones(Float32,N)
A_d = CUDA.ones(Float32,N)

@time wait(fookd(A_d, ndrange=N))
#@timev wait(fookh(A_h, ndrange=N))

# @kernel function foo(A)
#     I = @index(Global)
#     A[I]= I%UInt64
# end