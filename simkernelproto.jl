using KernelAbstractions
using KernelAbstractions.Extras
using CUDA
using CUDAKernels
using Random123
using WormlikeChain
using BenchmarkTools
using FastClosures

#BAOAB integration
@kernel function simulate!(positions, 
                          velocities,
                          externalforce!,
                          bondedforce,
                          @Const(beadparams),
                          invmass,
                          γ,
                          invβ,
                          VNbeads::Val{Nbeads},
                          step0::UInt64,
                          Vrngkey::Val{rngkey},
                          Δt,
                          VNdims::Val{Ndims},
                          VChainbounds::Val{Chainbounds},
                          VNsteps::Val{Nsteps}) where {Ndims, Nbeads, rngkey, Chainbounds, Nsteps}
    @uniform Nthreads = groupsize()[1]
    @uniform Nbeadsperthread= (Nbeads-1)÷Nthreads + 1
    @uniform Nrngsperbead= (Ndims-1)÷4 + 1
    tid = @index(Local, Linear)
    groupid = @index(Group, Linear)
    forces = @localmem Int64 (Nbeads, Ndims)
    #@inbounds sharedA[I] = 0
    #TODO set up force calc shared mem
    for step in 1:Nsteps
        #Reset forces
        #Kick
        @unroll for beadloopi in 1:Nbeadsperthread
            beadid= tid + (beadloopi-1)*Nthreads
            previd = WormlikeChain.prev_beadid(beadid,Chainbounds)
            nextid = WormlikeChain.next_beadid(beadid,Chainbounds) 
            #load bead1
            bead = ntuple(@closure(i -> @inbounds(positions[beadid,i])), Ndims)
            #load bead2
            #load bead3
            #load beadparams
            #calc bonded forces on bead1,2,3
            #force1, force2, force3= bondedforce(beadparam, bead1, bead2, bead3, step0+step)
            #inc force3
            @synchronize
            #inc force2
            @synchronize
            
            force = externalforce(beadparam, bead1)
            #externalforce!(forces,positions,beadparams,beadid)
            #apply total forces to each bead
            @unroll for i in 1:Ndims
                @inbounds velocities[beadid,i] += (invmass*Δt*(2f0^-32))*forces[beadid,i]
            end
        end
        #Drift with O step half way
        @unroll for beadloopi in 1:Nbeadsperthread
            mybeadid= tid + (beadloopi-1)*Nthreads
            randns= WormlikeChain.randn_32(rngkey, step0+step, (mybeadid<<32)%UInt64,VNdims)
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
pos_d = cu(pos_h)
vel_h = zeros(Float32,N,3)
vel_d = cu(vel_h) 
params_h = fill((k=4.0f0,),N)
params_d = cu(params_h)

function externalforce!(forces,positions,beadparams,id)
    @unroll for i in 1:3
        @inbounds f = -positions[id,i] * beadparams[id].k
        @inbounds forces[id,i] += unsafe_trunc(Int64,f*(2f0^32))
    end
end

function externalforce(beadparams,id)
    @unroll for i in 1:3
        @inbounds f = -positions[id,i] * beadparams[id].k
        @inbounds forces[id,i] += unsafe_trunc(Int64,f*(2f0^32))
    end
end



@btime wait(simkh(pos_h, 
                vel_h,
                externalforce!,
                params_h,
                1f0,
                1f0/10f0,
                1f0,
                Val(N),
                1%UInt64,
                Val(124432%UInt64),
                0.1f0,
                Val(3),
                Val((1,N+1)),
                Val(1000),
                 ndrange=N))

@btime wait(simkd(pos_d, 
                vel_d,
                externalforce!,
                params_d,
                1f0,
                1f0/10f0,
                1f0,
                Val(N),
                1%UInt64,
                Val(124432%UInt64),
                0.1f0,
                Val(3),
                Val((1,N+1)),
                Val(1000),
                 ndrange=N))