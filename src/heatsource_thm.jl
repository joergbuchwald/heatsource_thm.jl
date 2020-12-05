module heatsource_thm

#using Parameters
using SpecialFunctions

struct input_param
    E::Float32
    ν::Float32
    ρ_w::Float32
    ρ_s::Float32
    K_w::Float32
    K_s::Float32
    c_w::Float32
    c_s::Float32
    k::Float32
    μ::Float32
    a_s::Float32
    a_w::Float32
    T₀::Float32
    Q::Float32
    n::Float32
    aprime::Float32
    γ_w::Float32
    λ::Float32
    G::Float32
    K::Float32
    bprime::Float32
    m::Float32
    κ::Float32
    K_hydr::Float32
    a_u::Float32
    c::Float32
    X::Float32
    Y::Float32
    Z::Float32
end


function input_param(E, ν, rho_w, rho_s, K_w, K_s, c_w, c_s, k, mu, a_s, a_w, T0, Q, n)
    aprime = a_s
    γ_w = 9.81 * rho_w
    λ = E * ν / ((1+ν)*(1-2*ν))
    G = E / (2*(1+ν))
    K = n*K_w + (1-n)*K_s
    bprime = (λ+2*G/3)*aprime
    K_hydr = k*rho_w*9.81/mu
    m = n*rho_w*c_w+(1-n)*rho_s*c_s
    κ = K/m
    a_u = a_s*(1-n)+a_w*n
    c = K_hydr*(λ+2*G)/γ_w
    X = a_u*(λ+2*G)-bprime
    Y = 1/(λ+2*G) * (X/((1-c/κ)*a_u)+bprime/a_u)
    Z = 1/(λ+2*G) * (X/((1-c/κ)*a_u))
    input_param(E, ν, rho_w, rho_s, K_w, K_s, c_w, c_s, k, mu, a_s, a_w, T0, Q, n, aprime, γ_w, λ, G, K, bprime, m, κ, K_hydr, a_u, c, X, Y, Z)
end
# Test data:
param = input_param(5.0e9, 0.3, 999.1, 2290.0, 0.6, 1.838, 4280.0, 917.654, 2.0e-20, 1.0e-3, 1.5e-5, 4.0e-4, 273.15, 300.0, 0.16)

function R(x,y,z)
    return sqrt.(x.^2 .+ y.^2 .+ z.^2)
end

function f(κ, R, t)
    return erfc.(R ./ (2.0 .* sqrt.(κ.*t)))
end

function g(κ, R, t)
    return (κ .* t ./ R.^2 .+ (1.0 ./ 2. .- κ .* t ./ R.^2) .* erfc.(R ./ (2. .* sqrt.(κ .* t))) .- sqrt.(κ .* t ./ (pi .* R.^2)) .* exp.(-R.^2 ./ (4. .* κ .* t)))
end

function fstar(Y, Z, κ, c, R, t)
    return (Y*f(κ,R,t)-Z*f(c,R,t))
end

function gstar(Y, Z, κ, c, R, t)
    return (Y*g(κ,R,t)-Z*g(c,R,t))
end
function temperature(R, t, p::input_param)
    Q, K, κ, T₀ = p.Q, p.K, p.κ, p.T₀
    return (Q / (4*pi .* K .* R) * f(κ, R, t) .+ T₀)
end
function porepressure(R, t, p::input_param)
    X, Q, K, κ, c = p.X, p.Q, p.K, p.κ, p.c
    return (X / (1.0-c/κ) .* Q ./ (4.0 * pi * K .* R) .* (f(κ, R, t) .- f(c, R, t)))
end

function ux_i(R, x_i, t, p::input_param)
    a_u, Q, K, Y, Z, κ, c = p.a_u, p.Q, p.K, p.Y, p.Z, p.κ, p.c
    return a_u .* x_i .* Q/(4.0 * pi * K .* R) .* gstar(Y, Z, κ, c, R, t)
end

function dg_dR(κ, i, R, t)
    return ((2.0 .* i ./ R.^3) .* sqrt.(κ .* t ./ pi) .* exp.(-R .* R ./ (4.0 * κ .* t)) .+ (2.0 * i * κ .* t ./ R.^4) .* (f(κ, R, t) .- 1.0))
end
function dgstar_dR(Y, Z, κ, c, i, R, t) # Subscript R means derivative w.r.t R
    return (Y .* dg_dR(κ,i,R,t) .- Z .* dg_dR(c,i,R,t))
end

function sigma_N(x, y, z, ii, t, p::input_param) # N for normal components
    Q, a_u, K, G, Y, Z, κ, c, λ, bprime, T₀ = p.Q, p.a_u, p.K, p.G, p.Y, p.Z, p.κ, p.c, p.λ, p.bprime, p.T₀
    r = R(x,y,z)
    i = 0
    if ii == :xx
        i = x
    elseif ii == :yy
        i = y
    elseif ii == :zz
        i = z
    end
    return ((Q * a_u ./ (4.0 * pi * K .* r)) .* (2.0 * G*( gstar(Y, Z, κ, c, r, t) * (1 - i^2 ./ r^2) .+ i * dgstar_dR(Y, Z, κ, c, i,r,t))
                    .+ λ .* (x .* dgstar_dR(Y, Z, κ, c, x, r, t) .+ y .* dgstar_dR(Y, Z, κ, c, y, r, t) .+ z .* dgstar_dR(Y, Z, κ, c, z, r, t) .+
                    2.0 .* gstar(Y, Z, κ, c, r, t))) .- bprime .* (temperature(r,t,p) .- T₀))
end
function sigma_S(x, y, z, t, i, j, p::input_param) # S for shear components
    Q, a_u, K, G, Y, Z, κ, c, λ, bprime, T₀ = p.Q, p.a_u, p.K, p.G, p.Y, p.Z, p.κ, p.c, p.λ, p.bprime, p.T₀
    r = R(x,y,z)
    i = 0
    if i == :x
        i_ = x
    elseif i == :y
        i_ = y
    elseif i == :z
        i_ = z
    end
    if j == :x
        j_ = x
    elseif j == :y
        j_ = y
    elseif j == :z
        j_ = z
    end
    return ((Q * a_u ./ (4.0 * pi * K .* R)) .* (2.0 * G .* (i_ * dgstar_dR(Y, Z, κ, c, j_, R, t) / 2 .+ j_ .* dgstar_dR(Y, Z, κ, c, i_, R, t)
        / 2.0 .- i_ * j_ .* gstar(Y, Z, κ, c, R, t) / R.^2)))
end
end # module
