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

# ╔═╡ f655040a-3738-11eb-0bc3-c9ac7f870dff
struct point
	x::Float64
	y::Float64
	z::Float64
end

# ╔═╡ 3af4b608-3739-11eb-1ae6-5fabbf93fb13
pt = point(0.3, 0.5, 1.4)

# ╔═╡ 4c0dfe90-3739-11eb-3f94-4b491ab06c7e
t=10 .^ range( log10(0.1), log10(1e7), length = 101 ) 

# ╔═╡ 53fd229c-373c-11eb-157a-b3dbe5dce12e
python_sol = heatsource.ANASOL()

# ╔═╡ 77c1ea96-373a-11eb-19cd-136266a50500
p = heatsource_thm.param

# ╔═╡ 30376f8e-373a-11eb-250a-e17872572733
begin 
	plt_temp=plot()
	plot!(t,heatsource_thm.temperature(pt.x,pt.y,pt.z,t,p), label="julia")
	plot!(t,python_sol.temperature(pt.x,pt.y,pt.z,t), label="python")
	plt_temp
end

# ╔═╡ aad3ac82-373a-11eb-25e6-a1caf33133fa
begin
	plt_press = plot()
	plot!(t,heatsource_thm.porepressure(pt.x,pt.y,pt.z,t,p), label="julia")
	plot!(t,python_sol.porepressure(pt.x,pt.y,pt.z,t), label="python")
	plt_press
end

# ╔═╡ b9940756-373a-11eb-0dd9-6b571c2fd100
begin
	plt = plot()
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:x,p), label="u_x julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"x"), label="u_x python")
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:y,p), label="u_y julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"y"), label="u_x python")
	plot!(t,heatsource_thm.u_i(pt.x,pt.y,pt.z,t,:z,p), label="u_z julia")
	plot!(t,python_sol.u_i(pt.x,pt.y,pt.z,t,"z"), label="u_x python")
	plt
end

# ╔═╡ f8f6f868-373a-11eb-0b50-65b4fda0f6bb
begin
	plts1 = plot()
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:xx,p), label="σ_xx julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"xx"), label="σ_xx python")
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:yy,p), label="σ_yy julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"yy"), label="σ_yy python")
	plot!(t,heatsource_thm.sigma_ii(pt.x,pt.y,pt.z,t,:zz,p), label="σ_zz julia")
	plot!(t,python_sol.sigma_ii(pt.x,pt.y,pt.z,t,"zz"), label="σ_zz python")
	plts1
end

# ╔═╡ 6a131a22-373b-11eb-2a8a-f522221f1f60
begin
	plts2 = plot()
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:x,:y,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"x","y"), label="σ_xy python")
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:y,:z,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"y","z"), label="σ_yz python")
	plot!(t,heatsource_thm.sigma_ij(pt.x,pt.y,pt.z,t,:z,:x,p), label="σ_xy julia")
	plot!(t,python_sol.sigma_ij(pt.x,pt.y,pt.z,t,"z","x"), label="σ_zx python")
	plts2
end

# ╔═╡ Cell order:
# ╠═88a82b6a-3738-11eb-3799-1934287b0b2d
# ╠═4838192e-373a-11eb-0133-0d68a42eda21
# ╠═c32f07ec-373b-11eb-2657-0596fdf76e76
# ╠═30fe1b94-373c-11eb-3401-e74a8c57bec0
# ╠═fba355e4-3739-11eb-2455-939badcc0129
# ╠═f655040a-3738-11eb-0bc3-c9ac7f870dff
# ╠═3af4b608-3739-11eb-1ae6-5fabbf93fb13
# ╠═4c0dfe90-3739-11eb-3f94-4b491ab06c7e
# ╠═53fd229c-373c-11eb-157a-b3dbe5dce12e
# ╠═77c1ea96-373a-11eb-19cd-136266a50500
# ╠═30376f8e-373a-11eb-250a-e17872572733
# ╠═aad3ac82-373a-11eb-25e6-a1caf33133fa
# ╠═b9940756-373a-11eb-0dd9-6b571c2fd100
# ╠═f8f6f868-373a-11eb-0b50-65b4fda0f6bb
# ╠═6a131a22-373b-11eb-2a8a-f522221f1f60
