# use pkgs
using JuMP, Gurobi, LinearAlgebra, Test

# parameters
b0  = 0.9780549381510736
d0  = 0.2757743321951065

### BASE SOCP ###
# model
model = Model(Gurobi.Optimizer)

# variables
@variable(model, y[1:2])

# objective 
@objective(model, Min, y[1])

# constraint
@constraint(model, [y[1] - d0, y[2] - b0] in SecondOrderCone())

# optimize
## Note that resolving the problem, shows an zero objective in the solver print.
## This seems to be a Gurobi.jl problem
optimize!(model)
solution_summary(model)

# test the objective value
@test objective_value(model) ≈ 0.2757743321951065

### REFORMULATED SOCP ###

# n = cardinality of y
# k = number of constraints

# Reformulation of ||A(η)y + b(η)||₂ <= cᵀ(χ)y + d(χ) ∀ η ∈ 𝓩ˡ, χ ∈ 𝓩ʳ

# Right hand side pertubation set
# 𝓩ʳ = {χ : ∃ u : Pχ + Qu + p ∈ K}, where K is closed convex pointed cone or 
# polyhedral cone 
# => P, Q, p are needed for (a)
# Remark: if the cone K is a direct product of simpler cones K¹,..., Kˢ, it 
# takes the form 
# 𝓩ʳ = {χ : ∃ u¹,...,uᴿ : Pₛχ + Qᵣuʳ + pᵣ ∈ Kʳ, r = 1,...,R}
# Example: 𝓩 is an intersection of concetric co-axial box and ellipsiod,
# 𝓩 = {χ ∈ ℝᴸ : -1 ≤ χₗ ≤ 1, l ≤ L, √(∑ₗ₌₁₋ₗ (χₗ)²/(σₗ)²) ≤ Ω},
# where σₗ > 0 and Ω > 0 are given parameters, and becomes,
# 𝓩 = {χ ∈ ℝᴸ : P₁χ + p₁ ∈ K¹, P₂χ + p₂ ∈ K²},
# where P₁χ = [χ;0], p₁ = [zeros(L,1);1], K¹ = {(z,t) ∈ ℝᴸ × ℝ : t ≥ ||z||∞}, 
# whence its dual K¹* = {(z,t) ∈ ℝᴸ × ℝ : t ≥ ||z||₁}
# where P₂χ = [∑⁻¹χ;0] with ∑ = diagm(σ₁,..,σₗ), p₂ = [zeros(L,1);Ω] and K² is 
# the Lorentz cone of the dimension L+1, whence its dual K²* = K².

# Left hand side uncertainty set
# Zˡ ={η = [δA,δb] : |(δA)ᵢⱼ| ≤ δᵢⱼ, 1 ≤ i ≤ k, 1 ≤ j ≤ n, |(δb)ᵢ|, 1 ≤ i ≤ k},
# [A(ζ),b(ζ)] = [Aⁿ, bⁿ] + [δA, δb].
# => Aⁿ, bⁿ, δᵢⱼ, δᵢ are needed for (b)
# Example: if b ∈ [0.75,1.25], than bⁿ = 1.0 and corresponding δᵢ = 0.25

# Equivalent explicit system of conic quadratic and linear constraints
# (a.1) τ + pᵀv ≤ δ(y)
# (a.2) Pᵀv = σ(y)
# (a.3) Qᵀv = 0
# (a.4) v ∈ K* = {v : vᵀw ≥ 0, ∀ w ∈ K}
# (b.1) zᵢ ≥ |(Aⁿy + bⁿ)ᵢ| + δᵢ + ∑ⱼ₌₁₋ₙ |δᵢⱼyⱼ|, ∀ i ∈ 1,...,k 
# (b.2) ||z||₂ <= τ
# where K* is the dual cone of K, note that nonnegative orthants, Lorentz and 
# Semidefinite cones are self-dual, and thus their finite direct products, i.e.,
# canonical cones, are self-dual as wel. 

# additional variables
# v is a vector variable with the same cardinality of y
# z is a vector variable with the same cardinality of y
# τ is a scalar variable

# reformulation of |x| + |y| ≤ 1 => x + y ≤ 1, x - y ≤ 1, -x + y ≤ 1, -x - y ≤ 1

# model
rm = Model(Gurobi.Optimizer)

# variables
@variable(rm, y[1:2])
@variable(rm, v[1:2])
@variable(rm, z[1:2])

@variable(rm, τ)

# objective 
@objective(model, Min, y[1])

# constraints


# b2) ||z||₂ <= τ
@constraint(rm, [τ, [z...]] in SecondOrderCone())