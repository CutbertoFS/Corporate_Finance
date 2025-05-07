#= ################################################################################################## 

    Fin 971: Spring 2025 Corporate Finance
    Problem Set 7: DeMarzo and Fishman (2007, RFS)

    Last Edit:  May 6, 2025
    Authors:    Cutberto Frias Sarraf, Zachary Orlando

=# ##################################################################################################

using Parameters, Plots, Random, LinearAlgebra, Statistics, LaTeXStrings, Distributions, Serialization
using JuMP, Ipopt, Interpolations
using LaTeXStrings

#= ################################################################################################## 
    Parameters
=# ##################################################################################################

@with_kw struct Primitives

    β::Float64      = exp(-0.0953)                  # Investor's discount factor
    δ::Float64      = exp(-0.0998)                  # Agent's discount factor   
    L::Float64      = 75.0                          # Investor's termination value
    R::Float64      = 0.0                           # Agent's termination value 
    λ::Float64      = 1.0                           # Diverted fraction

    Y_L::Float64    = 0.0                           # Project cash flow Low value
    Y_H::Float64    = 20.0                          # Project cash flow High value
    π_H::Float64    = 0.5                           # Probability of a High cash flow value

    na::Int64       = 500

    μ::Float64      = π_H * Y_H + (1-π_H) * Y_L     # Expected cash flow

end 

# Initialize value function and policy functions
@with_kw mutable struct Results
    b::Vector{Float64}
    d_H::Vector{Float64}
    d_L::Vector{Float64}
    p_H ::Vector{Float64}
    p_L::Vector{Float64}
    a_H::Vector{Float64}
    a_L::Vector{Float64}
end

@with_kw struct OtherPrimitives
    a_max::Float64
    a_min::Float64   
    a_grid::Array{Float64,1}
end

# Function for initializing model primitives and results
function Initialize_Model()
    param = Primitives()
    @unpack_Primitives param

    a_max               = μ / (1 - β)
    a_min               = max(π_H * λ * (Y_H - Y_L) + R, δ^(-1) * R)
    # a_grid              = collect(range(a_min, length = na, stop = a_max))
    a_grid              = exp.(range(log(a_min), log(a_max), length=na))

    b                   = (μ - (β - δ) * δ^(-1) * R) / (1 - β) .- a_grid
    d_H                 = zeros(na)
    d_L                 = zeros(na)
    p_H                 = zeros(na)
    p_L                 = zeros(na)
    a_H                 = zeros(na)
    a_L                 = zeros(na)
    
    other_param         = OtherPrimitives(a_max, a_min, a_grid)
    results             = Results(b, d_H, d_L, p_H, p_L, a_H, a_L)

    return param, results, other_param
end

#= ################################################################################################## 
  
    Functions

=# ##################################################################################################

function Solver(a::Float64, old_b::Vector{Float64}, param::Primitives, results::Results, other_param::OtherPrimitives)
    @unpack_Primitives param
    @unpack_Results results
    @unpack_OtherPrimitives other_param

    Inter_func             = interpolate(old_b, BSpline(Cubic(Line(OnGrid()))))
    Interpolation_function = extrapolate(Inter_func, Line())

    # Define the model
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "print_level", 1)

    # Register the custom interpolation function with JuMP
    register(model, :custom_interp, 1, x -> Interpolation_function(clamp(x, a_min, a_max)); autodiff = true)

    # Decision variables
    @variable(model, 0 <= d_H <= a_max)
    @variable(model, 0 <= d_L <= a_max)
    @variable(model, 0 <= p_H <= 1)
    @variable(model, 0 <= p_L <= 1)
    @variable(model, a_min <= a_H <= a_max)
    @variable(model, a_min <= a_L <= a_max)
    set_start_value(d_H, 0.5 * a_max)
    set_start_value(d_L, 0.5 * a_max)
    set_start_value(p_H, 0.5)
    set_start_value(p_L, 0.5)
    set_start_value(a_H, 0.5 * (a_min + a_max))
    set_start_value(a_L, 0.5 * (a_min + a_max))

    # Interpolated continuation values
    @NLexpression(model, b_H, custom_interp(a_H))
    @NLexpression(model, b_L, custom_interp(a_L))

    # Promise-Keeping constraint
    @NLexpression(model, Payoff_Agent_H, (d_H + p_H * R) + δ * (1 - p_H) * a_H)
    @NLexpression(model, Payoff_Agent_L, (d_L + p_L * R) + δ * (1 - p_L) * a_L)
    @NLconstraint(model, a == π_H * Payoff_Agent_H + (1 - π_H) * Payoff_Agent_L)

    # Incentive Compatibility
    @NLexpression(model, IC_HH, d_H + (1 - p_H) * δ * a_H + p_H * R)
    @NLexpression(model, IC_HL, λ * (Y_H - Y_L) + d_L + (1 - p_L) * δ * a_L + p_L * R)
    @NLconstraint(model, IC_HH >= IC_HL)

    # Investor payoff
    @NLexpression(model, Payoff_Investor_H, (Y_H - d_H + p_H * L) + β * (1 - p_H) * b_H)
    @NLexpression(model, Payoff_Investor_L, (Y_L - d_L + p_L * L) + β * (1 - p_L) * b_L)
    @NLexpression(model, Objective, π_H * Payoff_Investor_H + (1 - π_H) * Payoff_Investor_L)
    @NLobjective(model, Max, Objective)

    optimize!(model)

    return (
        obj = objective_value(model),
        d_H = value(d_H),
        d_L = value(d_L),
        p_H = value(p_H),
        p_L = value(p_L),
        a_H = value(a_H),
        a_L = value(a_L),
        b_H = value(b_H),
        b_L = value(b_L)
    )
end


function Solve_Model(param::Primitives, results::Results, other_param::OtherPrimitives; tol::Float64 = 1e-3)
    @unpack_Primitives param
    @unpack_Results results
    @unpack_OtherPrimitives other_param

    next_b = zeros(na)
    old_b  = b
    error  = 10
    n      = 0

    while error > tol                                         # Begin iteration
        for i in 1:na
            next_b[i] = Solver(a_grid[i], old_b, param, results, other_param)[1]
            # println("Solver is in ", i, " iteration.")
        end 
        error = maximum(abs.(next_b .- old_b))                # Reset error level
        old_b = copy(next_b)                                  # Update value function
        n += 1
        println("Solve_Model is in ", n, " iteration with error ", error)
    end

    for i in 1:na
        res = Solver(a_grid[i], old_b, param, results, other_param)
        b[i], d_H[i], d_L[i], p_H[i], p_L[i], a_H[i], a_L[i] = res
    end

end

#= ################################################################################################## 
    Solving the Model
=# ##################################################################################################

param, results, other_param = Initialize_Model()
Solve_Model(param, results, other_param)
@unpack_Primitives param                                             
@unpack_Results results
@unpack_OtherPrimitives other_param

b_FB  = (μ - (β - δ) * δ^(-1) * R) / (1 - β) .- a_grid


# Value function for investor b(a)
plot(a_grid, b,
     label      = L"\mathsf{Asymmetric\ Info}",
     lw         = 2,                  
     color      = :black)
plot!(a_grid, b_FB,
     label      = L"\mathsf{First\ Best}",
     lw         = 2,                  
     color      = :blue)
xlabel!(L"\mathsf{Agents\ payoff\ a}")
ylabel!(L"\mathsf{Investors\ payoff\ b}")
plot!(legend = :topright,             
      xtickfont  = font(11),
      ytickfont  = font(11),
      guidefont  = font(12),
      legendfont = font(12),
      titlefont  = font(11),
      xlims  = (10, 62),
      ylims  = (30, 100),
      size       = (600, 600))   



plot(a_grid, b)
plot(a_grid, d_H)
plot(a_grid, d_L)
plot(a_grid, p_H)
plot(a_grid, p_L)
plot(a_grid, a_H)
plot(a_grid, a_L)