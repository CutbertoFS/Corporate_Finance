#= ################################################################################################## 

    Fin 971: Spring 2025 Corporate Finance
    Problem Set 7: DeMarzo and Fishman (2007, RFS)

    Last Edit:  May 6, 2025
    Authors:    Cutberto Frias Sarraf, Zachary Orlando

=# ##################################################################################################

using Parameters, Plots, Random, LinearAlgebra, Statistics, LaTeXStrings, Distributions, Serialization
using NLopt, Interpolations, Dierckx

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

    na::Int64       = 100

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

function Solver_NLopt(a::Float64, old_b::Vector{Float64}, param::Primitives, results::Results, other_param::OtherPrimitives)
    @unpack_Primitives param
    @unpack_Results results
    @unpack_OtherPrimitives other_param

    # Precompute interpolation function once
    Inter_func              = interpolate(old_b, BSpline(Cubic(Line(OnGrid()))))
    Interpolation_function  = extrapolate(Inter_func, Flat())

    function interp_b(a_H, a_L)
        b_H = Interpolation_function(clamp(a_H, a_min, a_max))
        b_L = Interpolation_function(clamp(a_L, a_min, a_max))
        return b_H, b_L
    end

    # 6 Decision Variables: [d_H, d_L, p_H, p_L, a_H, a_L]	
    Options = Opt(:LD_SLSQP, 6)
    lower_bounds!(Options, [0.0  , 0.0  , 0.0, 0.0, a_min, a_min])
    upper_bounds!(Options, [a_max, a_max, 1.0, 1.0, a_max, a_max])

    function objective(x::Vector, grad::Vector)
        d_H, d_L, p_H, p_L, a_H, a_L = x
        b_H, b_L                     = interp_b(a_H, a_L)
        PI_H                         = (Y_H - d_H + p_H * L) + β * (1 - p_H) * b_H
        PI_L                         = (Y_L - d_L + p_L * L) + β * (1 - p_L) * b_L
        return -(π_H * PI_H + (1 - π_H) * PI_L)
    end

    function constraints(result::Vector, x::Vector, grad::Matrix)
        d_H, d_L, p_H, p_L, a_H, a_L = x

        # 1. Promise-keeping
        PA_H        = (d_H + p_H * R) + δ * (1 - p_H) * a_H
        PA_L        = (d_L + p_L * R) + δ * (1 - p_L) * a_L
        result[1]   = π_H * PA_H + (1 - π_H) * PA_L - a

        # 2. Incentive compatibility
        IC_HH       = d_H + (1 - p_H) * δ * a_H + p_H * R
        IC_HL       = λ * (Y_H - Y_L) + d_L + (1 - p_L) * δ * a_L + p_L * R
        result[2]   = IC_HH - IC_HL
    end

    equality_constraint!(Options, constraints, [0.0, 0.0])

    xtol_rel!(Options, 1e-4)
    maxeval!(Options, 1000)

    x0 = [0.5 * a_max, 0.5 * a_max, 0.5, 0.5, 0.5 * (a_min + a_max), 0.5 * (a_min + a_max)]
    min_objective!(Options, objective)

    (value, solution, ret) = optimize(Options, x0)

    d_H, d_L, p_H, p_L, a_H, a_L = solution
    b_H, b_L                     = interp_b(a_H, a_L)

    return (
        obj = -value,
        d_H = d_H,
        d_L = d_L,
        p_H = p_H,
        p_L = p_L,
        a_H = a_H,
        a_L = a_L,
        b_H = b_H,
        b_L = b_L
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
            next_b[i] = Solver_NLopt(a_grid[i], old_b, param, results, other_param)[1]
            println("Solver is in ", i, " iteration.")
        end 
        error = maximum(abs.(next_b .- old_b))                # Reset error level
        old_b = copy(next_b)                               # Update value function
        n += 1
        println("Solve_Model is in ", n, " iteration with error ", error)
    end

    for i in 1:na
        res = Solver_NLopt(a_grid[i], old_b, param, results, other_param)
        b[i]    = res.obj
        d_H[i]  = res.d_H
        d_L[i]  = res.d_L
        p_H[i]  = res.p_H
        p_L[i]  = res.p_L
        a_H[i]  = res.a_H
        a_L[i]  = res.a_L
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
      xtickfont  = font(10),
      ytickfont  = font(10),
      guidefont  = font(13),
      legendfont = font(10),
      titlefont  = font(11),
      xlims  = (10, 62),
      ylims  = (30, 100),
      size       = (500, 500))   



plot(a_grid, b)
plot(a_grid, d_H)
plot(a_grid, d_L)
plot(a_grid, p_H)
plot(a_grid, p_L)
plot(a_grid, a_H)
plot!(a_grid, a_L)