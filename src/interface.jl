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
        end
        new(beaddef::BeadDefinition, perbead_params, global_params, init_pos, init_vel)
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
    SpecificForce(pe,force,params,interactions,VNdims,exprpe)
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
    chainbounds is a tuple of at least 2 Ints, 
        chainbounds[1] is 1
        chainbounds[end] is the total number of beads + 1
        the other element in chainbounds are the bead ids of the start of a new chain
    """
    chainbounds::Tuple{T,T,Vararg{T}} where T <: Integer
end

# function ChainSystem(chains, boundary_pe, added_global_params)
#     #symbolically diff potential energy to get force
#     Ndims= size(init_pos)[2]
#     @variables b[1:Ndims]
#     @variables bpars[1:Nbparams]
#     @variables gpars[1:Ngparams]
#     bsymparams= (; zip(perbead_param_keys,(bpars...,))...)
#     gsymparams= (; zip(global_param_keys,(gpars...,))...)
#     exprpe= chain_pe(p,b,n,bsymparams,gsymparams)
#     g= Symbolics.gradient(-exprpe, [p;b;n]; simplify=true)
#     @debug g
#     gexpr= build_function(g, [p;b;n], bpars, gpars ,expression=Val{true})
#     gfun= eval(gexpr[1])
#     function chain_force(p,b,n, bpars, gpars)
#         f= gfun([p;b;n], bpars, gpars)
#         pf = ntuple((i -> @inbounds(f[i])), VNdims)
#         bf = ntuple((i -> @inbounds(f[i+Ndims])), VNdims)
#         nf = ntuple((i -> @inbounds(f[i+2Ndims])), VNdims)
#         (SVector(pf),SVector(bf),SVector(nf))
#     end
#     BeadDefinition(chain_pe, chain_force, perbead_param_keys, global_param_keys, VNdims, exprpe)
# end