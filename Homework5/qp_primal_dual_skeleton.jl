################################################################################
## Code for solving equality-constrained Quadratic optimization problems with ##
## primal-dual IPM with the following steps:                                  ##
##                                                                            ##
## 1. Compute a stricly feasible primal solution                              ##
## 2. Compute a feasible dual solution by solving the analytic centering      ##
##    problem using Newton's method                                           ##
## 3. Starting from the analytic center, find the path to optimal solution    ##
##    using the primal-dual path-following IPM                                ##
##                                                                            ##
## NOTE: Your task is to complete the Newton direction update starting from   ##
##       line 178. Everything else is provided for in this code file.         ##
################################################################################
using JuMP
using Clp
using Plots
pyplot()
srand(1)

## Problem data in standard form
c = [-1.0,-1.0, 0.0,0.0,0.0,0.0,0.0]
b = [ 5.0,-1.0,-8.0,9.0,4.0]
A = [-1/3  1.0 1.0 0.0 0.0 0.0 0.0;
      1/5 -1.0 0.0 1.0 0.0 0.0 0.0;
     -8/3 -1.0 0.0 0.0 1.0 0.0 0.0;
      1/2  1.0 0.0 0.0 0.0 1.0 0.0;
      1.0 -1.0 0.0 0.0 0.0 0.0 1.0]
## Create random PSD matrix Q for the quadratic term
n = length(c)
Q = randn(n, n)                # Create random matrix
Q = (Q + Q')/2                 # Make A symmetric
if isposdef(A) == false        # Check if A is PD
    λmin = eigmin(Q)           # Minimum eigenvalue
    Q = Q + (abs(λmin) + 1)*I  # Add λmin + 1 to diagonal elements
end
Q[3:n,:] = 0.0                 # Set zero values for slack variables
Q[:,3:n] = 0.0

## Standard form:
# min      -x₁ - x₂ + ⁠(1/2)xᵀQx
# s.t.: -1/3x₁ + x₂ + x₃                 =  5
#        1/5x₁ - x₂     + x₄             = -1
#       -8/3x₁ - x₂         + x₅         = -8
#        1/2x₁ + x₂             + x₆     =  9
#           x₁ - x₂ +               + x₇ =  4
#           x₁, ... ,x₇ ≧  0


################################################################################
## Find an initial strictly feasible primal solution                          ##
################################################################################
function strict_primal(c::Vector{Float64}, A::Matrix{Float64}, b::Vector{Float64})

    (m,n) = size(A)

    lp = Model(solver = ClpSolver())

    @variable(lp, t >= 0)
    @variable(lp, x[1:n] >= 0)
    @constraint(lp, A*x .== b)
    @constraint(lp, [i = 1:n], t >= 1 - x[i])
    @objective(lp, Min, t)

    solve(lp)

    return (getvalue(t), getvalue(x))
end

################################################################################
## Compute the analytic center of the feasible region that gives both
## feasible primal and dual solutions. The analytic center is obtained
## by solving
##
##   minimize  f(x) = dot(c,x) + xᵀQx- μ*sum(log(x))
## subject to  Ax = b
##
## for large penalty value μ.
## We solve the problem with Newton's method using block elimination +
## exact line search.
################################################################################
function analytic_center(c, B, A, Q, f, x, N, μ; ϵ = 1e-6, λ = 1.0, β = 0.99)
    ## 1st and 2nd derivatives for line search
    D(θ, λ)  = ForwardDiff.derivative(θ, λ)
    D²(θ, λ) = ForwardDiff.derivative(λ -> ForwardDiff.derivative(θ, λ), λ)

    ## θ : line search function
    ## λᵢ: initial step size (e.g. 0)
    ## ϵ : tolerance
    ## β : step size reduction for feasibility
    ## d : direction
    ## x : primal solution
    ## Returns optimal step size
    function Newton_ls(θ, λᵢ, ϵ, β, d, x)
        ## Make sure not to have negative logarithm
        while minimum(x + λᵢ*d) <= 0
            λᵢ = β*λᵢ
        end
        while abs(D(θ, λᵢ)) > ϵ             # Check stopping condition
            λᵢ₊₁ = λᵢ - D(θ, λᵢ)/D²(θ, λᵢ)  # Update step size until optimal
            λᵢ   = λᵢ₊₁
        end
        return λᵢ
    end

    ## Solve analytic centering problem to get primal an dual feasible
    ## solutions using Newton's method with block elimination updates
    for i = 1:N-1
        H  =  μ*diagm(x[i,:].^-2) + Q       # Hessian
        ∇f =  c - μ.*x[i,:].^-1 + Q*x[i,:]  # Gradient of original objective
        v  = (A*inv(H)*A')\(A*inv(H)*∇f)    # Dual variable vector v
        u  =  μ.*x[i,:].^-1                 # Dual variable vector u
        if norm(∇f) < ϵ                     # Stopping condition #1
            return (x[1:i,:], v, u)         # Return primal + dual solutions
        end
        d  = -inv(H)*(∇f - A'*v)            # Newton direction
        ## Line search
        ls(λ) = f(x[i,:] + λ*d)
        λ = Newton_ls(ls, λ, ϵ, β, d, x[i,:])
        x[i+1,:] = x[i,:] + λ*d             # Move to a new point
        if i == N - 1                       # Stopping condition #2
            return (x, v, u)                # too few iterations
        end
    end
end

################################################################################
## Find initial feasible primal and dual solution                             ##
################################################################################
function initial_solution(μ::Float64, c::Vector{Float64},
    b::Vector{Float64}, A::Matrix{Float64}, Q::Matrix{Float64})

    n = length(c)  # Number of x and u variables
    N = 20         # Iterations for computing analytic center

    ## Logarithmic barrier function for computing analytic center
    f(x) = dot(c,x) + (1/2)*dot(x,Q*x) - μ*sum(log(x[i]) for i = 1:n)

    x      = zeros(N, n)             # Store values of x for plotting
    (t,x₀) = strict_primal(c, A, b)  # Get strictly feasible primal solution
    x[1,:] = x₀                      # Set initial solution

    (x₀,v₀,u₀) = analytic_center(c, b, A, Q, f, x, N, μ)

    ## Plot trajectory to analytic center
    x1 = linspace(0,15,100)
    x2 = linspace(0,15,100)
    plot(x1, 5 + (1/3)*x1,  color = :blue)
    plot!(x1, 1 + (1/5)*x1, color = :blue)
    plot!(x1, 8 - (8/3)*x1, color = :blue)
    plot!(x1, 9 - (1/2)*x1, color = :blue)
    plot!(x1, -4 + x1,      color = :blue,
        xaxis = ("x₁", (0,10)),
        yaxis = ("x₂", (0,10)),
        aspect_ratio = :equal,
        size = (800,800),
        legend = false)
    traj = x[1:end,1:2]
    plot!(traj[:,1], traj[:,2], marker = :o, title = "Centering step")
    gui()
    savefig("qp_analytic_center.pdf")

    return (v₀,u₀,x₀[end,:])
end

################################################################################
## Update the Newton direction                                                ##
################################################################################
function update_newton_dir(A::Matrix{Float64}, Q::Matrix{Float64}, u::Vector{Float64},
                           x::Vector{Float64}, e::Vector{Float64}, μ::Float64)

    (m,n) = size(A)
    ## Diagonalize u and x
    U = diagm(u)
    X = diagm(x)
    ## Update Newton direction for QP
    ## NOTE: Complete these based on the formulas you derived in 5.2 (c)
    dᵥ =
    dₓ =
    dᵤ =
    return (dᵥ,dᵤ,dₓ)
end

################################################################################
## Calculate step size that retains feasibility                             ##
################################################################################
function calculate_step_size(x::Vector{Float64}, d::Vector{Float64}, ϵ::Float64)
    n = length(d)
    α = 0.9999
    for i=1:n
        if d[i] < 0
            α = min(α, -x[i]/d[i]) # prevents variable becoming negative
        end
    end
    return round(α, Int(-log10(ϵ))) #rounding avoids numerical issues
end

################################################################################
## Solve the QP problem and its dual                                          ##
##                                                                            ##
## minimize    cᵀx + (1/2)xᵀQx        maximize    bᵀv - (1/2)xᵀQx             ##
## subject to  Ax = b                 subject to  Aᵀv + u - Qx = c            ##
##              x ≥ 0                                        u ≥ 0            ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
################################################################################
function primal_dual_ip(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64},
                        Q::Matrix{Float64}, β::Float64, ϵ::Float64, N::Int)

    n = length(c)    # Number of primal variables x and dual variables u
    m = length(b)    # Number of dual variables v
    x = zeros(N, n)  # Primal variable x values
    u = zeros(N, n)  # Dual variable u values
    v = zeros(N, m)  # Dual variable v values
    e = ones(n)      # Vector of ones
    α_p = 1          # Step size for variable x update
    α_d = 1          # Step size for variable u, and v update
    μ = 10000.0

    ## Find an initial solution
    (v₀,u₀,x₀) = initial_solution(μ, c, b, A, Q)

    v[1,:] = v₀
    u[1,:] = u₀
    x[1,:] = x₀

    ## Main loop
    for i = 1:N
        ## Stopping condition #1
        if n*μ < ϵ
            v = v[1:i,:]
            u = u[1:i,:]
            x = x[1:i,:]
            return (v,u,x)
        end
        ## Update dₓ, dᵥ, dᵤ
        (dᵥ,dᵤ,dₓ) = update_newton_dir(A,Q,u[i,:],x[i,:],e,μ)

        ## Calculate step size α primal and α dual
        α_p = calculate_step_size(x[i,:], dₓ, ϵ)
        α_d = calculate_step_size(u[i,:], dᵤ, ϵ)

        ## Update variables
        v[i+1,:] = v[i,:] + α_d*dᵥ
        u[i+1,:] = u[i,:] + α_d*dᵤ
        x[i+1,:] = x[i,:] + α_p*dₓ

        ## Stopping condition #2
        if i == N-1
            return (v,u,x)
        end
        ## Update μ
        μ = μ*β
    end
end

## Executions and plotting
N = 100     # Number of iterations
β = 0.2     # Reduction factor
ϵ = 1e-6    # Tolerance

## Solution time
@time (v,u,x) = primal_dual_ip(A, b, c, Q, β, ϵ, N)

## Plotting...
ng = 1000
x1 = linspace(0,15,ng)
x2 = linspace(0,15,ng)

plot(x1, 5 + (1/3)*x1,  color = :blue, reuse = false, legend = false)
plot!(x1, 1 + (1/5)*x1, color = :blue)
plot!(x1, 8 - (8/3)*x1, color = :blue)
plot!(x1, 9 - (1/2)*x1, color = :blue)
plot!(x1, -4 + x1, color = :blue,
    xaxis = ("x₁", (0,10)),
    yaxis = ("x₂", (0,10)),
    aspect_ratio = :equal,
    size = (800,800),
    legend = false,
    title = "Optimization step")

g(x) = dot(c[1:2],x) + (1/2)*dot(x,Q[1:2,1:2]*x)
z = [g([x1[i],x2[j]]) for j = 1:ng, i = 1:ng]

contour!(x1, x2, z,
        levels = [2, 9.27, 27, 64, 128, 256],
        clims = (0,300),
        clabels = true)

traj = x[1:end,1:2]
plot!(traj[:,1], traj[:,2], marker = :o)
gui()
savefig("qp_convergence.pdf")
