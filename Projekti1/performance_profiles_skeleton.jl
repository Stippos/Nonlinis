using Plots
pyplot()

## Create a file "solution_methods.jl" and make sure it's in the
## same folder as this file. Implement the following methods:
##
## - Newton's method
## - BFGS Quasi-Newton
## - Gradient descent
##
## Also, implement the following line search methods:
##
## - Newton's method (univariate)
## - Armijo's rule (give your own parameterization α and β)
##
## NOTE: Your task is to implement the methods and make the function
##       templates to return correct information in the following lines:
##
##       - 172, 174: Newton's method with Newton's line search
##       - 192, 194: Newton's method with Armijo's rule
##
##       - 212, 214: BFGS with Newton's line search
##       - 232, 234: BFGS with Armijo's rule
##
##       - 252, 254: Gradient descent with Newton's line search
##       - 272, 274: Gradient descent with Armijo's rule
##
## NOTE: DO NOT CHANGE THE FUNCTIONS OR OTHER PARAMETERS IN THIS FILE.
##
include("solution_methods.jl")

## Generate a random n x n symmetric positive definite matrix A and
## n x 1 vector b. Parameter δ controls the condition number of A
function generate_problem_data(n::Int, δ::Float64)
    A = randn(n,n)                 # Create random matrix
    A = (A + A')/2                 # Make A symmetric
    if isposdef(A) == false        # Check if A is PD
        λmin = eigmin(A)           # Minimum eigenvalue
        A = A + (abs(λmin) + δ)*I  # Add λmin + δ to diagonal elements
    end
    b = randn(n)                   # Create random vector b
    return (A,b)                   # Reulting matrix A is PD
end

## Generate k test instances of dimension n with different values δ
## to control the condition number of the matrices.
function generate_instances(k::Int, n::Int, δ::StepRangeLen)
    A = Dict{Int,Matrix{Float64}}()   # Store matrices A
    b = Dict{Int,Vector{Float64}}()   # Store vectors  b
    for i = 1:k
        (A[i], b[i]) = generate_problem_data(n, δ[i])
    end
    return (A, b)
end

srand(1)                 # Control randomness
k  = 100                 # Number of intances to generate
n  = 100                 # Dimension of n x n PD matrix A
δ1 = range(.05,.05, k)   # Moderate condition numbers for matrices A
δ2 = range(.01,.01, k)   # Larger condition numbers for matrices A

## Generate problem data with δ1 and δ2
(A1, b1) = generate_instances(k, n, δ1)
(A2, b2) = generate_instances(k, n, δ2)

## Function to minimize with two different data
f1(x,i) = (1/2)*dot(x, A1[i]*x) - dot(b1[i], x)
f2(x,i) = (1/2)*dot(x, A2[i]*x) - dot(b2[i], x)

## Optimal solution costs
fopt = zeros(k,2)
for i = 1:k
    x1 = A1[i]\b1[i]
    x2 = A2[i]\b2[i]
    fopt[i,1] = f1(x1,i)
    fopt[i,2] = f2(x2,i)
end

##### Preallocate arrays ########

# Solution costs
fval_newton_newton   = zeros(k,2)
fval_newton_armijo   = zeros(k,2)
fval_bfgs_newton     = zeros(k,2)
fval_bfgs_armijo     = zeros(k,2)
fval_gradient_newton = zeros(k,2)
fval_gradient_armijo = zeros(k,2)

# Solution times
time_newton_newton   = zeros(k,2)
time_newton_armijo   = zeros(k,2)
time_bfgs_newton     = zeros(k,2)
time_bfgs_armijo     = zeros(k,2)
time_gradient_newton = zeros(k,2)
time_gradient_armijo = zeros(k,2)

# Number of iterations
iter_newton_newton   = zeros(Int,k,2)
iter_newton_armijo   = zeros(Int,k,2)
iter_bfgs_newton     = zeros(Int,k,2)
iter_bfgs_armijo     = zeros(Int,k,2)
iter_gradient_newton = zeros(Int,k,2)
iter_gradient_armijo = zeros(Int,k,2)

# Solution status
stat_newton_newton   = fill(false,k,2)
stat_newton_armijo   = fill(false,k,2)
stat_bfgs_newton     = fill(false,k,2)
stat_bfgs_armijo     = fill(false,k,2)
stat_gradient_newton = fill(false,k,2)
stat_gradient_armijo = fill(false,k,2)

ns  = 6              # Number of solvers (methods) to compare
np  = k              # Number of problems to solve
tps = zeros(np,ns,2) # Computing times for each problem/method.
                     # tps[p,s,j] = Inf if method s fails to
                     # solve problem p.

N   = 1000          # Number of iterations
ϵ   = 1e-06          # Tolerance
x0  = ones(n)        # Starting point

## Go through all instances for both sets of data
for j = 1:2
    for i = 1:k
        println(i)
        g1(x) = f1(x,i)
        g2(x) = f2(x,i)
        ## NOTE: Implement the functions in the "solution_methods.jl"
        ##
        ## (fvalue, numiter) = Newton(g, x0, N, NEWTON, ϵ)    (line 172)
        ## (fvalue, numiter) = Newton(g, x0, N, ARMIJO, ϵ)    (line 192)
        ##
        ## (fvalue, numiter) = BFGS(g, x0, N, NEWTON, ϵ)      (line 212)
        ## (fvalue, numiter) = BFGS(g, x0, N, ARMIJO, ϵ)      (line 232)
        ##
        ## (fvalue, numiter) = Gradient(g, x0, N, NEWTON, ϵ)  (line 252)
        ## (fvalue, numiter) = Gradient(g, x0, N, ARMIJO, ϵ)  (line 272)
        ##
        ## Each function takes as argument:
        ##
        ##   g: function to minimize   (this is already provided)
        ##  x0: initial starting point (this is already provided)
        ##   N: number of iterations   (this is already provided)
        ##  LS: type of line search to use (NEWTON or ARMIJO)
        ##   ϵ: tolerance for stopping criterion (this is already provided)
        ##
        ## All of the arguments are already provided in this file
        ## except LS, which you can define as:
        ##
        ## @enum LS NEWTON ARMIJO
        ##
        ## in the julia file "solution_methods.jl". Then you can
        ## select which line search to use inside each function as:
        ##
        ## if LS == NEWTON
        ##     ## Code to call newton line search
        ## elseif LS == ARMIJO
        ##     ## Code to call Armijo line search
        ## end
        ##
        ## Each function returns
        ##
        ##   fvalue: objective value of the final solution
        ##  numiter: number of iterations
        ##
        ##  NOTE: Make sure to have a stopping criterion for
        ##        each method if the method has converged

        ## Newton + Newton's line search (LS)
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = Newton(g1, x0, N, NEWTON, ϵ)
        else
            (fvalue, numiter) = Newton(g2, x0, N, NEWTON, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_newton_newton[i,j] = fvalue       # Store objective value
        time_newton_newton[i,j] = soltime      # Store solution time
        iter_newton_newton[i,j] = numiter      # Store iteration count
        stat_newton_newton[i,j] = status       # Store solution status
        ## Set solution time accordingly
        if status == true
            tps[i,1,j] = soltime
        else
            tps[i,1,j] = Inf
        end

        ## Newton + Armijo
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = Newton(g1, x0, N, ARMIJO, ϵ)
        else
            (fvalue, numiter) = Newton(g2, x0, N, ARMIJO, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_newton_armijo[i,j] = fvalue       # Objective value
        time_newton_armijo[i,j] = soltime      # Solution time
        iter_newton_armijo[i,j] = numiter      # Iteration count
        stat_newton_armijo[i,j] = status       # Solution status
        ## Set solution time accordingly
        if status == true
            tps[i,2,j] = soltime
        else
            tps[i,2,j] = Inf
        end

        ## Gradient + Newton
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = BFGS(g1, x0, N, NEWTON, ϵ)
        else
            (fvalue, numiter) = BFGS(g2, x0, N, NEWTON, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_bfgs_newton[i,j] = fvalue         # Objective value
        time_bfgs_newton[i,j] = soltime        # Solution time
        iter_bfgs_newton[i,j] = numiter        # Iteration count
        stat_bfgs_newton[i,j] = status         # Solution status
        ## Set solution time accordingly
        if status == true
            tps[i,3,j] = soltime
        else
            tps[i,3,j] = Inf
        end

        ## Gradient + Armijo
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = BFGS(g1, x0, N, ARMIJO, ϵ)
        else
            (fvalue, numiter) = BFGS(g2, x0, N, ARMIJO, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_bfgs_armijo[i,j] = fvalue         # Objective value
        time_bfgs_armijo[i,j] = soltime        # Solution time
        iter_bfgs_armijo[i,j] = numiter        # Iteration count
        stat_bfgs_armijo[i,j] = status         # Solution status
        ## Set solution time accordingly
        if status == true
            tps[i,4,j] = soltime
        else
            tps[i,4,j] = Inf
        end

        ## BFGS + Newton
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = Gradient(g1, x0, N, NEWTON, ϵ)
        else
            (fvalue, numiter) = Gradient(g2, x0, N, NEWTON, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_gradient_newton[i,j] = fvalue     # Objective value
        time_gradient_newton[i,j] = soltime    # Solution time
        iter_gradient_newton[i,j] = numiter    # Iteration count
        stat_gradient_newton[i,j] = status     # Solution status
        ## Set solution time accordingly
        if status == true
            tps[i,5,j] = soltime
        else
            tps[i,5,j] = Inf
        end

        ## BFGS + Armijo
        starttime = time()                     # Start timer
        if j == 1
            (fvalue, numiter) = Gradient(g1, x0, N, ARMIJO, ϵ)
        else
            (fvalue, numiter) = Gradient(g2, x0, N, ARMIJO, ϵ)
        end
        soltime = time() - starttime           # Solution time
        status = abs(fvalue - fopt[i,j]) <= ϵ  # Check if solved or not
        fval_gradient_armijo[i,j] = fvalue     # Objective value
        time_gradient_armijo[i,j] = soltime    # Solution time
        iter_gradient_armijo[i,j] = numiter    # Iteration count
        stat_gradient_armijo[i,j] = status     # Solution status
        ## Set solution time accordingly
        if status == true
            tps[i,6,j] = soltime
        else
            tps[i,6,j] = Inf
        end
    end
end

for j = 1:2
    ###### Plot performance profiles ######
    tpsmin = minimum(tps[:,:,j], 2)  # Minimum time for each instance
    rps    = tps[:,:,j] ./ tpsmin    # Compute performance ratios
    τ      = sort(unique(rps))       # Sort the performance ratios in increasing order
    τ[end] == Inf && pop!(τ)         # Remove the Inf element if it exists

    ns = 6                           # Number of solvers
    np = k                           # Number of problems

    ρS = Dict()                      # Compute cumulative distribution functions
    for i = 1:ns                     # for performance ratios
        ρS[i] = [sum(rps[:,i] .<= τi) / np for τi in τ]
    end

    # Plot performance profiles
    labels = ["Newton (Exact)", "Newton (Armijo)", "BFGS (Exact)",
              "BFGS (Armijo)", "Gradient (Exact)", "Gradient (Armijo)"]
    styles = [:solid, :dash, :dot, :dashdot, :solid, :dash]
    fig = plot(xscale = :log2,   # OR :none, :ln, :log10
               xlabel = "τ",
               ylabel = "P(rₚₛ ≤ τ : 1 ≤ s ≤ nₛ)",
               title  = "Perfomance plot (condition $j)",
               size   = (1200,800),
               reuse  = false)
    for i = 1:ns
        plot!(τ, ρS[i], label = labels[i], seriestype = :steppre, linewidth = 2, line = styles[i])
    end
    ylims!(0, 1)
    xlims!(1, maximum(τ))
    savefig("performance_plot_$j.pdf") # Save figure as pdf
    gui()
end
