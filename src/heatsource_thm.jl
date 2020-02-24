module heatsource_thm

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
end

struct derived_param
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

p = input_param(2.8, 0.25, 999.1, 4280.0, 0.6, 1.838, 4280.0, 917.654, 2.0e-20, 1.0e-3, 1.5e-5, 4.0e-4, 273.15, 300, 0.16)

 
function determine_derivedconst(p::input_param)
    aprime = p.a_s
    γ_w = 9.81 * p.ρ_w 
    λ = p.E * p.ν / ((1+p.ν)*(1-2*p.ν))
    G = p.E / (2*(1+p.ν))
    K = p.n*p.K_w + (1-p.n)*p.K_s
    bprime = (λ+2*G/3)*aprime
    K_hydr = p.k*p.ρ_w*p.c_w+(1-p.n)*p.ρ_s*p.c_s
    m = p.n*p.ρ_w*p.c_w+(1-p.n)*p.ρ_s*p.c_s
    κ = K/m   
    a_u = p.a_s*(1-p.n)+p.a_w*p.n
    c = K_hydr*(λ+2*G)/γ_w
    X = a_u*(λ+2*G)-bprime
    Y = 1/(λ+2*G) * (X/((1-c/κ)*a_u)+bprime/a_u)
    Z = 1/(λ+2*G) * (X/((1-c/κ)*a_u))
    d = derived_param(aprime, γ_w, λ, G, K, bprime, m, κ, K_hydr, a_u, c, X, Y, Z)
    return d
end

d = determine_derivedconst(p)

function f(κ, R, t)
    return erfc(R/(2*sqrt(κ*t)))
end

function g(κ, R, t)
    return (κ*t/R^2+(1/2-κ*t/R**2)*erfc(R/(2*sqrt(κ*t)))-sqrt(κ*t/(pi*R^2))*exp(-R^2/(4*κ*t)))
end

function fstar(κ, R, t)
    return (Y*f(κ,R,t)-Z*f(c,R,t))
end

function gstart(κ, R, t)
    return (Y*g(κ,R,t)-Z*g(c,R,t))
end
end # module
