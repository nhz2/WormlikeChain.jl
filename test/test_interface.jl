@testset "BeadDefinition" begin
    param_keys= (:k,:kθ)
    function chain_pe(p,b,n,params)
        rp= (p .- b)
        rn= (n .- b)
        dp= sqrt(sum(rp.^2))
        dn= sqrt(sum(rn.^2))
        cosθ= sum(rp .* rn)/dp/dn
        params.k*1//2*((dn-1)^2) + params.kθ*(cosθ+1)
    end
    beaddef1= BeadDefinition(chain_pe, param_keys, Val(3))
    p= [0.0, 0.0, 0.0]
    b= [0.0, 0.0, 1.0]
    n= [0.0, 0.0, 2.0]
    params= (k=3.5, kθ=4.4)
    @test 0.0 == beaddef1.chain_pe(p,b,n,params)
    fp, fb, fn = beaddef1.chain_force(p,b,n,params)
    @test 0.0 == norm(fp)
    @test 0.0 == norm(fb)
    @test 0.0 == norm(fn)
end

@testset "Chain" begin
    param_keys= (:k,:kθ)
    function chain_pe(p,b,n,params)
        rp= (p .- b)
        rn= (n .- b)
        dp= sqrt(sum(rp.^2))
        dn= sqrt(sum(rn.^2))
        cosθ= sum(rp .* rn)/dp/dn
        params.k*1//2*((dn-1)^2) + params.kθ*(cosθ+1)
    end
    beaddef1= BeadDefinition(chain_pe, param_keys, Val(3))
    params= fill((k=3.5, kθ=4.4),10)
    params[10]= (k=0.0, kθ= 0.0)# turn off interation to split ends
    params[1]= (k=3.5, kθ= 0.0)
    pos= zeros(10,3)
    pos[:,3] .= range(0.0, step= 1.0, length=10)
    Chain(beaddef1,params, pos)
end