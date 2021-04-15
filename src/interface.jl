using Symbolics
using StaticArrays



"""
The potential and number of dimentions of a bead
"""
struct BeadDefinition{Ndims}
    chain_pe
    chain_force
    perbead_param_keys
    global_param_keys
    VNdims::Val{Ndims}
    chain_pe_expression
end

function BeadDefinition(chain_pe, perbead_param_keys::NTuple{Nbparams,Symbol}, global_param_keys::NTuple{Ngparams,Symbol} , VNdims::Val{Ndims}) where {Ndims, Nbparams, Ngparams}
    #symbolically diff potential energy to get force
    @variables p[1:Ndims]
    @variables b[1:Ndims]
    @variables n[1:Ndims]
    @variables bpars[1:Nbparams]
    @variables gpars[1:Ngparams]
    bsymparams= (; zip(perbead_param_keys,(bpars...,))...)
    gsymparams= (; zip(global_param_keys,(gpars...,))...)
    exprpe= chain_pe(p,b,n,bsymparams,gsymparams)
    g= Symbolics.gradient(-exprpe, [p;b;n]; simplify=true)
    @debug g
    gexpr= build_function(g, [p;b;n], bpars, gpars ,expression=Val{true})
    gfun= eval(gexpr[1])
    function chain_force(p,b,n, bpars, gpars)
        f= gfun([p;b;n], bpars, gpars)
        pf = ntuple((i -> @inbounds(f[i])), VNdims)
        bf = ntuple((i -> @inbounds(f[i+Ndims])), VNdims)
        nf = ntuple((i -> @inbounds(f[i+2Ndims])), VNdims)
        (SVector(pf),SVector(bf),SVector(nf))
    end
    BeadDefinition(chain_pe, chain_force, perbead_param_keys, global_param_keys, VNdims, exprpe)
end

struct Chain
    beaddef::BeadDefinition
    perbead_params
    global_params
    init_pos
    init_vel
    function Chain(beaddef::BeadDefinition, perbead_params, global_params, init_pos, init_vel)
        #error checking
        size(init_pos) == size(init_vel) || throw(DimensionMismatch("init_pos $(size(init_pos)) and init_vel $(size(init_vel)) don't match"))
        length(size(init_pos)) == 2 || throw(DimensionMismatch("init_pos $(length(size(init_pos))) must be 2"))
        Val(size(init_pos)[2]) == beaddef.VNdims || throw(DimensionMismatch("init_pos $(size(init_pos)) must match bead dimension $(beaddef.VNdims)"))
        Nbeads= size(init_pos)[1]
        Nbeads == length(perbead_params) || throw(DimensionMismatch("number of coords $(Nbeads) must match number of params $(length(perbead_params))"))
        #make sure forces are finite
        for beadi in 1:Nbeads
            p= init_pos[mod(beadi-2,Nbeads)+1,:]
            b= init_pos[beadi,:]
            n= init_pos[mod(beadi,Nbeads)+1,:]
            pf, bf, nf = beaddef.chain_force(p,b,n,perbead_params[beadi], global_params)
            e= beaddef.chain_pe(p,b,n,perbead_params[beadi], global_params)
            all(isfinite.([pf;bf;nf;e])) || @warn "some forces or energies are not finite on Bead number $(beadi)"
            all(abs.([pf;bf;nf;e]) .< 2.0^30) || @warn "some forces or energies are over 2^30 on Bead number $(beadi)"
        end
        new(beaddef::BeadDefinition, copy(perbead_params), global_params, copy(init_pos), copy(init_vel))
    end
end

"""
zero velocities by default
"""
function Chain(beaddef::BeadDefinition, perbead_params, global_params, init_pos)
    Chain(beaddef::BeadDefinition, perbead_params, global_params, init_pos, zero(init_pos))
end

"""
A force that acts on a specific list of beads
"""
struct SpecificForce{Ndims,Nways}
    pe
    force
    params
    interactions::AbstractArray{NTuple{Nways,NTuple{2,T}},1} where T<:Integer
    VNdims::Val{Ndims}
    pe_expression
end

"""
Create a specific force
pe 

"""
function SpecificForce(pe, params, 
    interactions::AbstractArray{NTuple{Nways,NTuple{2,T}},1} where T<:Integer,
    VNdims::Val{Ndims}) where {Ndims, Nways}
    #symbolically diff potential energy to get force
    a = Num(SymbolicUtils.Sym{Real}(:a))
    @variables r[1:Nways,1:Ndims]
    @variables pars[1:length(params)]
    @variables time
    symparams= (; zip(keys(params),(pars...,))...)
    exprpe= pe(ntuple(i->r[i,:],Nways), time ,symparams) + 0time#This 0time is to force exprpe into a Num type
    rssymb= vcat(ntuple(i->r[i,:],Nways)...)
    g= Symbolics.gradient(-exprpe, rssymb; simplify=true)
    @debug g
    gexpr= build_function(g, rssymb, time, pars, expression=Val{true})
    gfun= eval(gexpr[1])
    function force(rs, time, pars)
        f= gfun(vcat(rs...), time, pars)
        ntuple(j->SVector(ntuple(i->f[i+Ndims*(j-1)],Ndims)),Nways)
    end
    SpecificForce(pe,force,params,copy(interactions),VNdims,exprpe)
end


"""
A collection of Chains and SpecificForces
"""
struct ChainSystem
    starttime
    specificforces
    beaddef::BeadDefinition
    perbead_params
    global_params
    init_pos
    init_vel
    """
    chainbounds is a tuple of at least 1 Int, 
        chainbounds[1] is 1
        chainbounds[end] is the total number of beads + 1
        the other element in chainbounds are the bead ids of the start of a new chain
    """
    chainbounds::Tuple{T,Vararg{T}} where T <: Integer
end

"""
Create an empty System
"""
function ChainSystem(starttime, beaddef::BeadDefinition, global_params)
    ChainSystem(starttime,(),beaddef,nothing,global_params,nothing,nothing,(1,))
end




"""
return a new ChainSystem with added chain 
"""
function append(s::ChainSystem, c::Chain)
    #error checking
    s.beaddef == c.beaddef || error("system and chain must have matching beaddef")
    s.global_params == c.global_params || error("system and chain must have matching global_params")
    Nbeads= size(c.init_pos)[1]
    if isnothing(s.init_pos)
        init_pos= c.init_pos
        init_vel= c.init_vel
        perbead_params= c.perbead_params
    else
        init_pos= [s.init_pos; c.init_pos]
        init_vel= [s.init_vel; c.init_vel]
        perbead_params= [s.perbead_params; c.perbead_params]
    end
    chainbounds= (s.chainbounds...,s.chainbounds[end]+Nbeads)
    news= ChainSystem(s.starttime,s.specificforces,s.beaddef,perbead_params,s.global_params,init_pos, init_vel,chainbounds)
    #check energy and forces for issues
    force_pe(news,news.init_pos,news.starttime)
    news
end

"""
return a new ChainSystem with added chain 
"""
function append(s::ChainSystem, f::SpecificForce)
    #error checking
    s.beaddef.VNdims == f.VNdims || error("system and force must have matching beaddef")
    specificforces= (s.specificforces..., f)
    news= ChainSystem(s.starttime,specificforces,s.beaddef,s.perbead_params,s.global_params,s.init_pos, s.init_vel,s.chainbounds)
    #check energy and forces for issues
    force_pe(news,news.init_pos,news.starttime)
    news
end


"""
Return force, pe of the beads
"""
function force_pe(s::ChainSystem,pos,time)
    @assert size(s.init_pos) == size(pos)
    f= zero(pos)
    pe= 0.0
    for sf in s.specificforces
        for interaction in sf.interactions
            rs= map(ids->pos[tobeadid(ids...,s.chainbounds),:], interaction)
            fs= sf.force(rs,time,sf.params)
            spe= sf.pe(rs,time,sf.params)
            all(isfinite.(vcat(fs..., spe))) || @warn "some forces or energies are not finite on specific force $(sf), interaction $(interaction)"
            all(abs.(vcat(fs..., spe)) .< 2.0^30) || @warn "some forces or energies are over 2^30 on specific force $(sf), interaction $(interaction))"
            for i in 1:length(interaction)
                f[tobeadid(interaction[i]...,s.chainbounds),:] .+= fs[i]
            end
            pe += spe
        end
    end
    Nbeads= size(pos)[1]
    #chain forces
    for beadid in 1:Nbeads
        previd = prev_beadid(beadid,s.chainbounds)
        nextid = next_beadid(beadid,s.chainbounds) 
        p= pos[previd,:]
        b= pos[beadid,:]
        n= pos[nextid,:]
        pf, bf, nf = s.beaddef.chain_force(p,b,n,s.perbead_params[beadid], s.global_params)
        e= s.beaddef.chain_pe(p,b,n,s.perbead_params[beadid], s.global_params)
        all(isfinite.([pf;bf;nf;e])) || @warn "chain forces or energies are not finite on Bead number $(beadid)"
        all(abs.([pf;bf;nf;e]) .< 2.0^30) || @warn "chain forces or energies are over 2^30 on Bead number $(beadid)"
        pe +=e
        f[previd,:] .+= pf
        f[beadid,:] .+= bf
        f[nextid,:] .+= nf
    end
    f, pe
end
