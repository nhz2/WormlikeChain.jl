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
                          numbeads::Val{Nbeads},
                          stepnum::UInt64,
                          rngkey::UInt64,
                          Δt,
                          beadsize::Val{Ndims},
                          chainstarts::Val{Chstarts},
                          steps::Val{Nsteps}) where {Ndims, Nbeads, Chstarts, Nsteps}
    @uniform Nthreads = groupsize()[1]
    @uniform Nbeadsperthread= (Nbeads-1)÷Nthreads + 1
    @uniform Nrngsperbead= (Ndims-1)÷4 + 1
    tid = @index(Local, Linear)
    groupid = @index(Group, Linear)

    #sharedA = @localmem Float32 N
    #@inbounds sharedA[I] = 0
    #TODO set up force calc shared mem
    #for step in 1:Nsteps
        #Kick
        #Drift with O step half way
        @unroll for beadloopi in 1:Nbeadsperthread
            mybeadid= tid + (beadloopi-1)*Nthreads + (groupid-1)*Nbeads
            @inbounds mypos= tuple((positions[mybeadid,i] for i in 1:Ndims)...)
            @inbounds myvel= tuple((velocities[mybeadid,i] for i in 1:Ndims)...)
            mypos = mypos .+ (Δt/2) .* myvel
            # velocities[mybeadid,:] .*= exp(-Δt*γ)
            # randns= WormlikeChain.randn_32(rngkey, stepnum, (mybeadid<<32)%UInt64,Ndims)
            # velocities[mybeadid,:] .+= (√(1-exp(-2γ*Δt))*√(invβ*invmass)) .* randns
            #mypos = mypos .+ (Δt/2) .* myvel
            #for i in 1:Ndims; positions[mybeadid,i]= mypos[i]; end
            #for i in 1:Ndims; velocities[mybeadid,i]= myvel[i]; end
        end
        #@synchronize
    #end
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
#                 Val(1),
#                  ndrange=N))

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
                Val(1),
                 ndrange=N))
