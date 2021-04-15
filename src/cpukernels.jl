

"""
Update positions and velocities with a simple leapfrog langevin dynamics BBAOA integrator 
"""
function refcpukernel!(pos, 
    vel,
    specificforces,
    specificinteractions,
    specificparams,
    chain_force,
    perbead_params,
    global_params,
    invmass,
    γ,
    invβ,
    VNbeads::Val{Nbeads},
    step0::UInt64,
    time0,
    Vrngkey::Val{rngkey},
    Δt,
    VNdims::Val{Ndims},
    VChainbounds::Val{Chainbounds},
    VNsteps::Val{Nsteps}) where {Ndims, Nbeads, rngkey, Chainbounds, Nsteps}
    force= zero(vel)
    for step in 1:Nsteps
        time= (step-1)*Δt + time0
        force .= 0
        #add specific forces
        for (sf, si, sp) in zip(specificforces,specificinteractions,specificparams)
            for interaction in si
                rs= map(beadid->pos[beadid,:], interaction)
                fs= sf(rs,time,sp)
                for i in 1:length(interaction)
                    force[interaction[i],:] .+= fs[i]
                end
            end
        end
        #add chain forces and do kick/ drift BAOAB stuff
        for beadid in 1:Nbeads
            previd = prev_beadid(beadid,Chainbounds)
            nextid = next_beadid(beadid,Chainbounds) 
            p= pos[previd,:]
            b= pos[beadid,:]
            n= pos[nextid,:]
            pf, bf, nf = chain_force(p,b,n,perbead_params[beadid], global_params)
            force[previd,:] .+= pf
            force[beadid,:] .+= bf
            force[nextid,:] .+= nf
            vel[beadid,:] .+= (invmass*Δt) .* force[beadid,:]
            pos[beadid,:] .+= (Δt/2) .* vel[beadid,:]
            vel[beadid,:] .*= exp(-Δt*γ)
            randns= randn_32(rngkey, step0+step, (beadid<<32)%UInt64,VNdims)
            vel[beadid,:] .+= (√(1-exp(-2γ*Δt))*√(invβ*invmass)) .* randns
            pos[beadid,:] .+= (Δt/2) .* vel[beadid,:]
        end
    end
end