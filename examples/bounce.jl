### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9285d6f0-9e71-11eb-0411-1928431042b7
begin
	import Pkg
	Pkg.activate(mktempdir())
end

# ╔═╡ f455bcc9-02bb-4d20-961d-4da33624a575
begin
	Pkg.add("Plots")
	Pkg.add("PlutoUI")
	Pkg.add("StaticArrays")
	Pkg.add("BenchmarkTools")
	using Plots
	using PlutoUI
end

# ╔═╡ d6c5e8aa-2b42-4ce1-9b42-e6c54bdcbf91
using WormlikeChain

# ╔═╡ 5a21fa3e-7053-4b85-860d-d011cb1f88cd
using StaticArrays

# ╔═╡ 18d60643-0f9e-430b-bff1-86aa191f8eae
using BenchmarkTools

# ╔═╡ f4577ff0-f627-4f54-87d3-27f4518e51db
md"""
## Defining a bead
A bead is defined by a potential energy function on it's neighbors and its dimension
"""

# ╔═╡ 362e0a98-ea1a-49b9-a3b4-a50dcdf6d10a
perbead_param_keys= (:k,:kθ)

# ╔═╡ cb945668-43cf-48aa-893c-2b968af8549e
function chain_pe(p,b,n,perbead_params,global_params)
	rp= (p .- b)
	rn= (n .- b)
	dp= sqrt(sum(rp.^2))
	dn= sqrt(sum(rn.^2))
	cosθ= sum(rp .* rn)/dp/dn
	perbead_params.k*1//2*((dn-1)^2) + perbead_params.kθ*(cosθ+1)
end

# ╔═╡ 4e5f3d5e-0cab-4fea-905d-4be4acef6526
const Ndims=2

# ╔═╡ cca6e1b6-84a2-4a1e-ab1c-a31d4b8ecb74
beaddef= BeadDefinition(chain_pe,  perbead_param_keys, (), Val(Ndims))

# ╔═╡ 3a27b649-d3f5-47b7-860d-fb67f1ac9f1f
md"""
BeadDefinition automatically differentiates `chain_pe` to get a `chain_force` 
`chain_force` can then be compiled into optimized assembly

"""

# ╔═╡ 4242953f-0f86-47a1-beab-bd2583a3dbd5
md"""
## Defining a chain
A chain is just a circle of beads.
"""

# ╔═╡ 09e78f9a-20f4-4ba6-8cb5-018897e96be1
const Nbeads=50

# ╔═╡ 8a6b422a-f288-4cbf-b8f0-ca9d41d8282a
begin
	r= (Nbeads/2π)
	θs= range(0.0, step= 2π/Nbeads, length=Nbeads)
	xs= (Nbeads/2π) .* cos.(θs)
	ys= (Nbeads/2π) .* sin.(θs)
	pos= [xs ys]
	# pos= zeros(Nbeads,Ndims)
	# pos[:,2] .= range(0.0, step= 1.0, length=Nbeads)
	vel= zero(pos)# .+ [0.0 -0.1]
	scatter(pos[:,1],pos[:,2])
	#beadparams= fill((k=0.0, kθ=0.0),Nbeads)
	beadparams= fill((k=3.5, kθ=4.4),Nbeads)
	# beadparams[Nbeads]= (k=0.0, kθ= 0.0)# turn off interation to split ends
	# beadparams[1]= (k=3.5, kθ= 0.0)
end

# ╔═╡ e12c77f3-870c-40ea-aec8-ae01750e26ee
chain1= Chain(beaddef,beadparams, (), pos,vel)

# ╔═╡ 440db016-35da-4e7e-9255-af9c42995a13
md"""
## Defining an additional potential energy term
"""

# ╔═╡ 78bfb5e5-9ec3-4c9f-9a43-329bf2087fb8
wallparams= (k= 1.0, L= 3*Nbeads/π)

# ╔═╡ 21a861fc-06d6-45bd-9a18-5f3573582693
function harmonicwallpe(pos,time,params)
	pe= 0
	for x in (pos[1])
		pe+= 1//2*params.k*(max(x-params.L/2,0))^2
		pe+= 1//2*params.k*(min(x+params.L/2,0))^2
	end
	pe
end

# ╔═╡ 06a89ca3-2768-4551-85bd-b24160b9b0e8
interactions = [((1,i),) for i in 1:Nbeads]

# ╔═╡ 3736c303-1e9c-44f2-b92a-09e632c19bfd
walldef= SpecificForce(harmonicwallpe, wallparams, interactions, Val(Ndims))

# ╔═╡ ed818b30-5e04-472b-b54f-6d2d5e77280f
md"""
## Putting it all together in a ChainSystem
First create an empty system, then add in chains and extra forces
"""

# ╔═╡ 7b971a08-e603-47e9-bf73-1c53b2f15d38
empty_system= ChainSystem(0.0,beaddef,())

# ╔═╡ c84c7c2b-1791-44ed-861d-b7076f3311b2
chain1system= append(empty_system,chain1)

# ╔═╡ fa9b09f7-6ab7-44a1-9e88-9ec1e108510f
system= append(chain1system,walldef)

# ╔═╡ b2f82994-caed-42dd-8fb8-3ad965446e35
sum(force_pe(system, system.init_pos,0.0)[1][:,2])

# ╔═╡ 4cd004e3-beeb-47a9-94a2-a4738c637ecc
randn(10,2)

# ╔═╡ eb54332d-3f8e-4884-92b1-949267412239
md"""
## Compiling the kernel and running a sim
Kernel get recompiled and speciallized for the new types created
"""

# ╔═╡ a88ade2a-66a6-443d-91a3-d48526905c86
totalframes= 10000

# ╔═╡ c4a9e6bb-befe-4403-9e97-46a8b43fe970
begin
	Δt= 0.01
	γ= 0.0
	invβ= 0.01
	trajectory= zeros(Nbeads,Ndims,totalframes)
	pes= zeros(totalframes)
	kes= zeros(totalframes)
	simpos= copy(system.init_pos)
	simvel= copy(system.init_vel)
	simvel[1,:] = [0.1 0.1]
	simvel .+= [0 1]
	for step in 1:totalframes
		trajectory[:,:,step] .= simpos
		f, pe= force_pe(system,simpos,0.0)
		pes[step] = pe
		kes[step] = 1//2*sum(simvel .^ 2)
		simvel .+= Δt.*f
		simpos .+= Δt.*simvel
		simvel .*= exp(-Δt*γ)
		simvel .+= (√(1-exp(-2γ*Δt))*√(invβ)) .* randn(Nbeads,Ndims)
		simpos .+= Δt.*simvel
		f, pe= force_pe(system,simpos,0.0)
		simvel .+= Δt.*f
	end
	nothing
end

# ╔═╡ 587def7f-a48b-4140-8571-bcf9151e0409
@bind seeframe Slider(1:totalframes)

# ╔═╡ 036a7fd3-11a0-4aae-a5d5-b819c1bf114b
scatter(trajectory[:,1,seeframe],trajectory[:,2,seeframe]; xlims = (-wallparams.L/2, wallparams.L/2), ylims = (-wallparams.L/2,wallparams.L/2))

# ╔═╡ 39fe110d-20b0-4362-b4d7-1f51a9fbf128
seeframe

# ╔═╡ e3db67ad-11f6-42f1-aaad-d81f3a772397
begin
	plot(kes)
	plot!(pes)
	plot!(kes+pes)
end

# ╔═╡ c3c2b1b4-4e73-4b68-8f2d-39fe1db6f550
plotly()

# ╔═╡ Cell order:
# ╟─f4577ff0-f627-4f54-87d3-27f4518e51db
# ╠═362e0a98-ea1a-49b9-a3b4-a50dcdf6d10a
# ╠═cb945668-43cf-48aa-893c-2b968af8549e
# ╠═4e5f3d5e-0cab-4fea-905d-4be4acef6526
# ╠═cca6e1b6-84a2-4a1e-ab1c-a31d4b8ecb74
# ╟─3a27b649-d3f5-47b7-860d-fb67f1ac9f1f
# ╟─4242953f-0f86-47a1-beab-bd2583a3dbd5
# ╠═09e78f9a-20f4-4ba6-8cb5-018897e96be1
# ╠═8a6b422a-f288-4cbf-b8f0-ca9d41d8282a
# ╠═e12c77f3-870c-40ea-aec8-ae01750e26ee
# ╟─440db016-35da-4e7e-9255-af9c42995a13
# ╠═78bfb5e5-9ec3-4c9f-9a43-329bf2087fb8
# ╠═21a861fc-06d6-45bd-9a18-5f3573582693
# ╠═06a89ca3-2768-4551-85bd-b24160b9b0e8
# ╠═3736c303-1e9c-44f2-b92a-09e632c19bfd
# ╟─ed818b30-5e04-472b-b54f-6d2d5e77280f
# ╠═7b971a08-e603-47e9-bf73-1c53b2f15d38
# ╠═c84c7c2b-1791-44ed-861d-b7076f3311b2
# ╠═fa9b09f7-6ab7-44a1-9e88-9ec1e108510f
# ╠═b2f82994-caed-42dd-8fb8-3ad965446e35
# ╠═4cd004e3-beeb-47a9-94a2-a4738c637ecc
# ╟─eb54332d-3f8e-4884-92b1-949267412239
# ╠═c4a9e6bb-befe-4403-9e97-46a8b43fe970
# ╠═a88ade2a-66a6-443d-91a3-d48526905c86
# ╠═036a7fd3-11a0-4aae-a5d5-b819c1bf114b
# ╠═39fe110d-20b0-4362-b4d7-1f51a9fbf128
# ╠═587def7f-a48b-4140-8571-bcf9151e0409
# ╠═e3db67ad-11f6-42f1-aaad-d81f3a772397
# ╠═9285d6f0-9e71-11eb-0411-1928431042b7
# ╠═f455bcc9-02bb-4d20-961d-4da33624a575
# ╠═c3c2b1b4-4e73-4b68-8f2d-39fe1db6f550
# ╠═d6c5e8aa-2b42-4ce1-9b42-e6c54bdcbf91
# ╠═5a21fa3e-7053-4b85-860d-d011cb1f88cd
# ╠═18d60643-0f9e-430b-bff1-86aa191f8eae
