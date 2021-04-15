@testset "BeadDefinition" begin
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
    p= [0.0, 0.0, 0.0]
    b= [0.0, 0.0, 1.0]
    n= [0.0, 0.0, 2.0]
    params= (k=3.5, kθ=4.4)
    globalparams= (d₀=1.0, )
    @test 0.0 == beaddef1.chain_pe(p,b,n,params, globalparams)
    fp, fb, fn = beaddef1.chain_force(p,b,n,params, globalparams)
    @test 0.0 == norm(fp)
    @test 0.0 == norm(fb)
    @test 0.0 == norm(fn)
end

@testset "Chain" begin
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
    params= fill((k=3.5, kθ=4.4),10)
    params[10]= (k=0.0, kθ= 0.0)# turn off interation to split ends
    params[1]= (k=3.5, kθ= 0.0)
    pos= zeros(10,3)
    pos[:,3] .= range(0.0, step= 1.0, length=10)
    globalparams= (d₀=1.0, )
    Chain(beaddef1,params, globalparams, pos)
end


@testset "SpecificForce" begin
    #set parameters for the force 
    params= (k=10.0,)
    #define potential energy, note x[1] is the 1st bead
    pe(x,time,params)= 1//2 * params.k * sum((x[1]).^2)
    #define the list of beads affected by this force
    #The inner tuple is the chain id and bead id in that chain
    #in this case, The force is external, 
    #so only one chain/bead id is listed
    interactions = [((1,1),)
                    ((1,2),)]
    #create the force Val(3) specifies 3d beads
    fdef= SpecificForce(pe, params, interactions, Val(3))
    @test [-10,-10,-10] == fdef.force(([1.0,1.0,1.0],), 0.0, params)[1]

    #now an example of a two bead time dependent interaction
    params2= (k₁=10.0, k₂=2.0)
    #note that pe2 is a regular function so it has local scope
    pe2(x,time,params)= 1//2 * sum((params.k₁, params.k₂, params.k₂*time) .* ((x[1] .- x[2]).^2))
    #define the list of beads affected by this force
    #The inner tuple is the chain id and bead id in that chain.
    #In this case, the force is on two beads, 
    #so only two chain/bead ids are listed per row.
    #In this example chain 1 bead 1 is bonded to chain 2 bead 1, and
    #chain 1 bead 2 is bonded to chain 1 bead 1
    interactions2 = [((1,1),(2,1))
                    ((1,2),(1,1))]
    #create the force Val(3) specifies 3d beads
    fdef2= SpecificForce(pe2, params2, interactions2, Val(3))
    #second argument is time
    f1,f2 = fdef2.force(([1.0,1.0,1.0],[0.0,0.0,0.0]), 0.0, params2)
    @test f1 == [-10,-2,0]
    @test f2 == .-f1
    #check force at time 1.0
    f1,f2 = fdef2.force(([1.0,1.0,1.0],[0.0,0.0,0.0]), 1.0, params2)
    @test f1 == [-10,-2,-2]
    @test f2 == .-f1

    #now test that SpecificForce will catch basic pe errors
    pe_bad(x,time,params)= 1//2 * params.k * (x[1][4])^2
    @test_throws BoundsError SpecificForce(pe_bad, params, interactions, Val(3))
end


@testset "ChainSystem" begin
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

end
