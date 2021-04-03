using KernelAbstractions
using KernelAbstractions.Extras
using CUDA
using CUDAKernels
using Random123
using WormlikeChain
using BenchmarkTools

#BAOAB integration
@kernel function simulate!(positions, 
                          velocities,
                          invmass,
                          γ,
                          invβ,
                          VNbeads::Val{Nbeads},
                          stepnum::UInt64,
                          Vrngkey::Val{rngkey},
                          Δt,
                          VNdims::Val{Ndims},
                          chainstarts::Val{Chstarts},
                          VNsteps::Val{Nsteps}) where {Ndims, Nbeads, rngkey, Chstarts, Nsteps}
    @uniform Nthreads = groupsize()[1]
    @uniform Nbeadsperthread= (Nbeads-1)÷Nthreads + 1
    @uniform Nrngsperbead= (Ndims-1)÷4 + 1
    tid = @index(Local, Linear)
    groupid = @index(Group, Linear)

    #sharedA = @localmem Float32 N
    #@inbounds sharedA[I] = 0
    #TODO set up force calc shared mem
    for step in 1:Nsteps
        #Kick
        #Drift with O step half way
        @unroll for beadloopi in 1:Nbeadsperthread
            mybeadid= tid + (beadloopi-1)*Nthreads + (groupid-1)*Nbeads
            randns= WormlikeChain.randn_32(rngkey, stepnum+step-1, (mybeadid<<32)%UInt64,VNdims)
            @unroll for i in 1:Ndims
                @inbounds positions[mybeadid,i] += (Δt/2) * velocities[mybeadid,i]
                @inbounds velocities[mybeadid,i] *= exp(-Δt*γ)
                @inbounds velocities[mybeadid,i] += (√(1-exp(-2γ*Δt))*√(invβ*invmass)) * randns[i]
                @inbounds positions[mybeadid,i] += (Δt/2) * velocities[mybeadid,i]
            end
        end
        @synchronize
    end
end

N=1024

simkh = simulate!(CPU(), N)
simkd = simulate!(CUDADevice(), N)

pos_h = zeros(Float32,N,3)
pos_d = CUDA.zeros(Float32,N,3)
vel_h = zeros(Float32,N,3)
vel_d = CUDA.zeros(Float32,N,3)

# @time wait(simkh(pos_h, 
#                 vel_h,
#                 1f0,
#                 1f0/10f0,
#                 1f0,
#                 Val(N),
#                 1%UInt64,
#                 124432%UInt64,
#                 1f0,
#                 Val(3),
#                 Val((1,)),
#                 Val(1000),
#                  ndrange=N))

@btime wait(simkd(pos_d, 
                vel_d,
                1f0,
                1f0/10f0,
                1f0,
                Val(N),
                1%UInt64,
                Val(124432%UInt64),
                1f0,
                Val(3),
                Val((1,)),
                Val(1000),
                 ndrange=N))
