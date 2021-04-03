using KernelAbstractions
using KernelAbstractions.Extras
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
                          VNbeads::Val{Nbeads},
                          stepnum::UInt64,
                          rngkey::UInt64,
                          Δt,
                          VNdims::Val{Ndims},
                          chainstarts::Val{Chstarts},
                          VNsteps::Val{Nsteps}) where {Ndims, Nbeads, Chstarts, Nsteps}
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
            randns= WormlikeChain.randn_32(rngkey, stepnum, (mybeadid<<32)%UInt64,VNdims)
            for i in 1:Ndims
                positions[mybeadid,i] += (Δt/2) * velocities[mybeadid,i]
                velocities[mybeadid,i] *= exp(-Δt*γ)
                velocities[mybeadid,i] += (√(1-exp(-2γ*Δt))*√(invβ*invmass)) * randns[i]
                positions[mybeadid,i] += (Δt/2) * velocities[mybeadid,i]
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

@time wait(simkh(pos_h, 
                vel_h,
                1f0,
                1f0/10f0,
                1f0,
                Val(N),
                1%UInt64,
                124432%UInt64,
                1f0,
                Val(3),
                Val((1,)),
                Val(1000),
                 ndrange=N))

@time wait(simkd(pos_d, 
                vel_d,
                1f0,
                1f0/10f0,
                1f0,
                Val(N),
                1%UInt64,
                124432%UInt64,
                1f0,
                Val(3),
                Val((1,)),
                Val(1000),
                 ndrange=N))
