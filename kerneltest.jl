using KernelAbstractions
using KernelAbstractions.Extras
using CUDA
using CUDAKernels
using Random123
using WormlikeChain
using StaticArrays

const N= 256

function getrow(A,row)
    ntuple(i -> A[row,i], Val(3))
end


@kernel function foo(A,B)
    I = @index(Global)
    myA= getrow(A,I)
    myB= getrow(B,I)
    myA= myA .+ myB
    @inbounds A[I,1] +=  myA[1]
    
end

fookh = foo(CPU(), 8)
fookd = foo(CUDADevice(), N)


A_h = ones(Float32,N,3)
A_d = CUDA.ones(Float32,N,3)
B_h = ones(Float32,N,3)
B_d = CUDA.ones(Float32,N,3)

#@time wait(fookd(A_d, B_d, ndrange=N))
#@time wait(fookh(A_h, B_h, ndrange=N))
show(A_d)