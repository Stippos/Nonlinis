using JuMP
using Ipopt
using ForwardDiff
using Plots
pyplot()

## Gradient and Hessian
∇(f,x)  = ForwardDiff.gradient(f,x)
∇²(f,x) = ForwardDiff.hessian(f,x)

## Problem data
n = 2  # Number of variables
m = 4  # Number of constraints (and dual variables)

## Objective function, its gradient and Hessian
f(x)   = 2x[1]^2 + 2x[2]^2 - 2(x[1]*x[2]) - 4x[1] - 6x[2]
∇f(x)  = ∇(f,x)
∇²f(x) = ∇²(f,x)

## Constraint functions
g₁(x) = x[1]^2 - x[2]     # ≦ 0
g₂(x) = x[1] + 5x[2] - 5  # ≦ 0
g₃(x) = -x[1]             # ≦ 0
g₄(x) = -x[2]             # ≦ 0

## Constraint function gradients and Hessians
∇g₁(x)  = ∇(g₁,x)
∇²g₁(x) = ∇²(g₁,x)       # Only g₁ has nonzero Hessian
∇g₂(x)  = ∇(g₂,x)
∇g₃(x)  = ∇(g₃,x)
∇g₄(x)  = ∇(g₄,x)

## Lagrangian function of the poblem
L(x,u) = f(x) + u[1]*g₁(x) + u[2]*g₂(x) + u[3]*g₃(x) + u[4]*g₄(x)
∇²L(x,u) = ∇²f(x) + u[1]*∇²g₁(x)  # g₂, g₃, and g₄ vanish as they are linear.

## Initialize data
k = 1             # Iteration count
N = 10            # Max iterations
ϵ = 1e-6          # Tolerance
xᵏ = [0.5; 0.5]   # Initial primal solutionm
uᵏ = zeros(4)     # initial dual solution
x = zeros(2, N)   # Save trajectory of iteratation
x[:,k] = xᵏ       # Set initial value

## Main loop
for k = 1:N-1
    xᵏ = x[:,k]
    ## Compute ∇f, ∇g, and ∇L²
    ∇fᵏ  = ∇f(xᵏ)
    ∇gᵏ  = [∇g₁(xᵏ), ∇g₂(xᵏ), ∇g₃(xᵏ), ∇g₄(xᵏ)]
    ∇L²ᵏ = ∇²L(xᵏ,uᵏ)
    gᵏ   = [g₁(xᵏ), g₂(xᵏ), g₃(xᵏ), g₄(xᵏ)]

    ## Projected lagrangian subproblem (Direction search)
    QP = Model(solver=IpoptSolver())
    @variable(QP, d[1:n])
    @objective(QP, Min, vecdot(∇fᵏ,d) + 0.5*(d'*∇L²ᵏ*d))
    @constraint(QP, LinearIneq[i=1:m], vecdot(∇gᵏ[i],d) + gᵏ[i] <= 0)

    solve(QP)

    dᵏ = getvalue(d)          # Obtain new direction
    x[:,k+1] = xᵏ + dᵏ        # Update primal solution
    uᵏ = getdual(LinearIneq)  # Obtain optimal dual solution

    ## Check stopping condition
    if norm(dᵏ) < ϵ
        x = x[:,1:k+1]
        break
    end
end

## Plotting
n  = 1000
x1 = linspace(-1,2,n)
x2 = linspace(-1,2,n)
z  = [f([x1[i],x2[j]]) for j = 1:n, i = 1:n]

## Contours of the objective
contour(x1, x2, z,
        levels  = [-11, -9, -7, -5, -3, -1],
        clims = (-20,5),
        clabels = true,
        legend = false)
## Plot the feasible region
plot!(x1, x1.^2, fill = (10,0.2), color = 2)
plot!(x1, 1 - (1/5)*x1, fill = (0,0.2), color = 3,
      xaxis = ("x₁", (0,2)),
      yaxis = ("x₂", (0,2)),
      aspect_ratio = :equal,
      size = (800,800),
      legend = false)
## Plot the trejectory of iterations
plot!(x[1,:], x[2,:], marker = :o, color = 1)

## First order approximation at initial point of constraint 1
plot!(x1, -0.25 .+ x1, color = 3, line = :dash)

savefig("SQP_example.pdf")
gui()
