using KernelAbstractions
using CUDA
using CUDAKernels
using Random123
using WormlikeChain

#BAOAB integration
@kernel function simulate!(positions, 
                          velocities,
                          invmass,
                          γ,
                          invβ,
                          numbeads::Val{Nbeads},
                          stepnum::UInt64,
                          rngkey::UInt64,
                          Δt
                          ;beadsize::Val{Ndims}=Val(3),
                          chainstarts::Val{Chstarts}=Val((1,)),
                          steps::Val{Nsteps}=Val(1)) where {Ndims, Nbeads, Chstarts, Nsteps}
    tid = @index(Local, Linear)
    groupid = @index(group, Linear)
    @uniform utid= tid%UInt64
    Nthreads = @uniform groupsize()[1]
    @uniform Nbeadsperthread= (Nbeads-1)÷Nthreads + 1
    @uniform Nrngsperbead= (Ndims-1)÷4 + 1
    #sharedA = @localmem Float32 N
    #@inbounds sharedA[I] = 0
    #TODO set up force calc shared mem
    for step in 1:Nsteps
        #Kick
        @synchronize
        #Drift with O step half way
        for i in 1:Nbeadsperthread
            mybeadid= tid + (i-1)*Nthreads
            positions[mybeadid,:] .+= (Δt*1//2) .* velocities[mybeadid,:]
            for rngstep in 1:Nrngsperbead
                randns= WormlikeChain.randn_4x32(rngkey, stepnum, mybeadid*Nrngsperbead-rngstep)
            velocities[mybeadid,:] .*= exp(-Δt*γ)
            velocities[mybeadid,:] .+= (√(1-exp(-2γ*Δt))*√(invβ*invmass)) .* (randns[1:end-1])
            positions[mybeadid,:] .+= (Δt*1//2) .* velocities[mybeadid,:]

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
