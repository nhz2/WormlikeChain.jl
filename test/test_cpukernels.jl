@testset "example kernel run" begin
    #make system
    perbead_param_keys= (:k,:kθ)
    global_param_keys= (:d₀,)
    function chain_pe(p,b,n,perbead_params,global_params)
        rp= (p .- b)
        rn= (n .- b)
        dp= sqrt(sum(rp.^2))
        dn= sqrt(sum(rn.^2))
        cosθ= sum(rp .* rn)/dp/dn
        perbead_params.k*1//2*((dn-global_params.d₀)^2) + perbead_params.kθ*(cosθ+1)
    end
    beaddef1= BeadDefinition(chain_pe, perbead_param_keys, global_param_keys, Val(3))
    globalparams= (d₀=1.0, )
    system= ChainSystem(0.0,beaddef1,globalparams)
    #create some chain
    params= fill((k=3.5, kθ=4.4),10)
    params[10]= (k=0.0, kθ= 0.0)# turn off interation to split ends
    params[1]= (k=3.5, kθ= 0.0)
    pos= zeros(10,3)
    pos[:,3] .= range(0.0, step= 1.0, length=10)
    chain1= Chain(beaddef1,params, globalparams, pos)
    system= append(system,chain1)
    #create chain2 translated [1 0 0] from chain1
    chain2= Chain(beaddef1,params, globalparams, pos .+ [1 0 0])
    system= append(system,chain2)
    #specific forces
    params= (k=10.0,)
    pe(x,time,params)= 1//2 * params.k * sum((x[1]).^2)
    interactions = [((1,1),)
                    ((1,2),)]
    fdef= SpecificForce(pe, params, interactions, Val(3))
    system= append(system,fdef)

    params2= (k₁=10.0, k₂=2.0)
    pe2(x,time,params)= 1//2 * sum((params.k₁, params.k₂, params.k₂*time) .* ((x[1] .- x[2]).^2))
    interactions2 = [((1,1),(2,1))
                    ((1,2),(1,1))]
    fdef2= SpecificForce(pe2, params2, interactions2, Val(3))
    system= append(system,fdef2)


    simpos= copy(system.init_pos)
    simvel= copy(system.init_vel)
    invmass= 1.0
    γ= 0.0
    invβ= 1.0
    VNbeads= Val(size(simpos)[1])
    VNdims= Val(size(simpos)[2])
    step0= UInt64(0)
    time0= 0.0
    Vrngkey= Val(UInt64(123))
    Δt= 0.01
    VChainbounds= Val(system.chainbounds)
    VNsteps= Val(2)
    specificinteractions= map(sf -> sf.interactions, system.specificforces)
    specificinteractions= map(arr -> map(beads2ids -> map(beads2id-> WormlikeChain.tobeadid(beads2id...,system.chainbounds),beads2ids), arr),specificinteractions)


    refcpukernel!(simpos, 
                simvel,
                map(sf->sf.force, system.specificforces),
                specificinteractions,
                map(sf->sf.params, system.specificforces),
                system.beaddef.chain_force,
                copy(system.perbead_params),
                system.global_params,
                invmass,
                γ,
                invβ,
                VNbeads,
                step0,
                time0,
                Vrngkey,
                Δt,
                VNdims,
                VChainbounds,
                VNsteps)
end