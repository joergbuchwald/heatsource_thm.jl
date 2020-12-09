### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 4838192e-373a-11eb-0133-0d68a42eda21
using Plots

# ╔═╡ c32f07ec-373b-11eb-2657-0596fdf76e76
using PyCall

# ╔═╡ fba355e4-3739-11eb-2455-939badcc0129
using heatsource_thm

# ╔═╡ 88a82b6a-3738-11eb-3799-1934287b0b2d
begin
    import DarkMode
    DarkMode.enable()
end

# ╔═╡ 30fe1b94-373c-11eb-3401-e74a8c57bec0
@pyimport heatsource

# ╔═╡ 588e72e0-3742-11eb-1ac0-6b670d8fed99
@pyimport ogs

# ╔═╡ 168247ec-394e-11eb-1b0d-a3c0f544e4f7
@pyimport vtuIO

# ╔═╡ f655040a-3738-11eb-0bc3-c9ac7f870dff
struct point
	x::Float64
	y::Float64
	z::Float64
end

# ╔═╡ 3af4b608-3739-11eb-1ae6-5fabbf93fb13
pt = point(0.3, 1.4, 0.0)

# ╔═╡ 4c0dfe90-3739-11eb-3f94-4b491ab06c7e
t=10 .^ range( log10(0.1), log10(1e7), length = 101 ) 

# ╔═╡ 53fd229c-373c-11eb-157a-b3dbe5dce12e
python_sol = heatsource.ANASOL()

# ╔═╡ 77c1ea96-373a-11eb-19cd-136266a50500
p = heatsource_thm.param

# ╔═╡ a4399ee2-3952-11eb-3520-cd3615442718
ogs6_data = vtuIO.PVDIO("/home/buchwalj/github/heatsource_thm.jl/", "blubb.pvd", dim=2)

# ╔═╡ 06323a1e-3953-11eb-2728-532c26e99ce1
pts = Dict("pt0" => (pt.x, pt.y, pt.z))

# ╔═╡ 2dd6d6ce-3953-11eb-0cfa-bddc2c177205
ogs6_press = ogs6_data.readTimeSeries("pressure_interpolated", pts=pts)

# ╔═╡ aad3ac82-373a-11eb-25e6-a1caf33133fa
begin
	plt_press = plot()
	plot!(t,heatsource_thm.porepressure(pt.x,pt.y,pt.z,t,p), label="julia")
	plot!(t,python_sol.porepressure(pt.x,pt.y,pt.z,t), label="python")
	plot!(ogs6_data.timesteps, ogs6_press["pt0"], label="ogs6")
	plt_press
end

# ╔═╡ 74a0b9c6-3953-11eb-1afa-31722dac8762
ogs6_temp = ogs6_data.readTimeSeries("temperature_interpolated", pts=pts)

# ╔═╡ 30376f8e-373a-11eb-250a-e17872572733
begin 
	plt_temp=plot()
	plot!(t,heatsource_thm.temperature(pt.x,pt.y,pt.z,t,p), label="julia")
	plot!(t,python_sol.temperature(pt.x,pt.y,pt.z,t), label="python")
	plot!(ogs6_data.timesteps, ogs6_temp["pt0"], label="ogs6")
	plt_temp
end

# ╔═╡ 748680ea-3978-11eb-15a0-f785c6d849e0
ogs6_displ = ogs6_data.readTimeSeries("displacement", pts=pts)

# ╔═╡ a284c36c-3a57-11eb-26e1-bf4d287be275
begin
	ogs6_ux=[]
	ogs6_uy=[]
end

# ╔═╡ b9940756-373a-11eb-0dd9-6b571c2fd100
begin
	plt = plot()
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:x,p), label="u_x julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"x"), label="u_x python")
	plot!(ogs6_data.timesteps,ogs6_ux, label="u_x ogs6")
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:y,p), label="u_y julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"y"), label="u_y python")
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:z,p), label="u_z julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"z"), label="u_z python")
	plot!(ogs6_data.timesteps,ogs6_uy, label="u_y ogs6")
	plt
end

# ╔═╡ 7658abde-3a57-11eb-3c26-0370131ab504
for i in 1:length(ogs6_displ["pt0"])
	push!(ogs6_ux,ogs6_displ["pt0"][i][1])
	push!(ogs6_uy,ogs6_displ["pt0"][i][2])
end

# ╔═╡ 91178dc2-3a58-11eb-0f64-73b38273e740


# ╔═╡ 902bab8e-3978-11eb-2d2f-95a57671124e
ogs6_sigma = ogs6_data.readTimeSeries("sigma", pts=pts)

# ╔═╡ 9ceff3f4-3978-11eb-33a2-4fde48d1acbe
begin
	ogs6_σxx=[]
	ogs6_σyy=[]
	ogs6_σzz=[]
	ogs6_σxy=[]
end

# ╔═╡ f8f6f868-373a-11eb-0b50-65b4fda0f6bb
begin
	plts1 = plot()
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:xx,p), label="σ_xx julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"xx"), label="σ_xx python")
	plot!(ogs6_data.timesteps,ogs6_σxx, label="σ_xx ogs6")
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:yy,p), label="σ_yy julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"yy"), label="σ_yy python")
	plot!(ogs6_data.timesteps,ogs6_σyy, label="σ_yy ogs6")
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:zz,p), label="σ_zz julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"zz"), label="σ_zz python")
	plot!(ogs6_data.timesteps,ogs6_σzz, label="σ_zz ogs6")
	plts1
end

# ╔═╡ 6a131a22-373b-11eb-2a8a-f522221f1f60
begin
	plts2 = plot()
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:x,:y,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"x","y"), label="σ_xy python")
	plot!(ogs6_data.timesteps,ogs6_σxy, label="σ_xy ogs6")
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:y,:z,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"y","z"), label="σ_yz python")
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:z,:x,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"z","x"), label="σ_zx python")
	plts2
end

# ╔═╡ 2e453874-3a59-11eb-1486-593d6afbb10b
for i in 1:length(ogs6_sigma["pt0"])
	push!(ogs6_σxx,ogs6_sigma["pt0"][i][1])
	push!(ogs6_σyy,ogs6_sigma["pt0"][i][2])
	push!(ogs6_σzz,ogs6_sigma["pt0"][i][3])
	push!(ogs6_σxy,ogs6_sigma["pt0"][i][4])
end

# ╔═╡ 66ad5fd0-3742-11eb-2acd-43b5fbfee758
model=ogs.OGS(INPUT_FILE="nummodel.prj", PROJECT_FILE="nummodel.prj")

# ╔═╡ a584f0ce-3742-11eb-2d17-73a152ff28ff
function writeinput()
	model.replaceParameter(name="E", value=p.E)
	model.replaceParameter(name="nu", value=p.ν)
	model.replaceParameter(name="T0", value=p.T₀)
	model.replaceParameter(name="temperature_ic", value=p.T₀)
	model.replaceParameter(name="temperature_bc_left", value=p.T₀)
	model.replaceParameter(name="temperature_source_term", value=p.Q/2)
	model.replacePhaseProperty(mediumid=0, phase="AqueousLiquid", name="specific_heat_capacity", value=p.c_w)
	model.replacePhaseProperty(mediumid=0, phase="AqueousLiquid", name="thermal_conductivity", value=p.K_w)
	model.replacePhaseProperty(mediumid=0, phase="AqueousLiquid", name="thermal_expansivity", value=p.a_w)
	model.replacePhaseProperty(mediumid=0, phase="AqueousLiquid", name="viscosity", value=p.μ)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="permeability", value=p.k)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="porosity", value=p.n)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="density", value=p.ρ_s)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="thermal_conductivity", value=p.K_s)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="specific_heat_capacity", value=p.c_s)
	model.replacePhaseProperty(mediumid=0, phase="Solid", name="thermal_expansivity", value=p.a_s/3)
	model.writeInput()
end

# ╔═╡ 527b13da-3743-11eb-03ad-a5f79e7f2eea
writeinput()

# ╔═╡ 6009b1c8-374d-11eb-1645-45351e981a8c
#model.runModel()

# ╔═╡ Cell order:
# ╠═88a82b6a-3738-11eb-3799-1934287b0b2d
# ╠═4838192e-373a-11eb-0133-0d68a42eda21
# ╠═c32f07ec-373b-11eb-2657-0596fdf76e76
# ╠═30fe1b94-373c-11eb-3401-e74a8c57bec0
# ╠═588e72e0-3742-11eb-1ac0-6b670d8fed99
# ╠═168247ec-394e-11eb-1b0d-a3c0f544e4f7
# ╠═fba355e4-3739-11eb-2455-939badcc0129
# ╠═f655040a-3738-11eb-0bc3-c9ac7f870dff
# ╠═3af4b608-3739-11eb-1ae6-5fabbf93fb13
# ╠═4c0dfe90-3739-11eb-3f94-4b491ab06c7e
# ╠═53fd229c-373c-11eb-157a-b3dbe5dce12e
# ╠═77c1ea96-373a-11eb-19cd-136266a50500
# ╠═a4399ee2-3952-11eb-3520-cd3615442718
# ╠═30376f8e-373a-11eb-250a-e17872572733
# ╠═aad3ac82-373a-11eb-25e6-a1caf33133fa
# ╠═b9940756-373a-11eb-0dd9-6b571c2fd100
# ╠═f8f6f868-373a-11eb-0b50-65b4fda0f6bb
# ╠═6a131a22-373b-11eb-2a8a-f522221f1f60
# ╠═06323a1e-3953-11eb-2728-532c26e99ce1
# ╠═2dd6d6ce-3953-11eb-0cfa-bddc2c177205
# ╠═74a0b9c6-3953-11eb-1afa-31722dac8762
# ╠═748680ea-3978-11eb-15a0-f785c6d849e0
# ╠═a284c36c-3a57-11eb-26e1-bf4d287be275
# ╠═7658abde-3a57-11eb-3c26-0370131ab504
# ╠═91178dc2-3a58-11eb-0f64-73b38273e740
# ╠═902bab8e-3978-11eb-2d2f-95a57671124e
# ╠═9ceff3f4-3978-11eb-33a2-4fde48d1acbe
# ╠═2e453874-3a59-11eb-1486-593d6afbb10b
# ╠═66ad5fd0-3742-11eb-2acd-43b5fbfee758
# ╠═a584f0ce-3742-11eb-2d17-73a152ff28ff
# ╠═527b13da-3743-11eb-03ad-a5f79e7f2eea
# ╠═6009b1c8-374d-11eb-1645-45351e981a8c
