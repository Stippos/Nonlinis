using JuMP
using Clp   #### Use Clp to compare LP solution oF ADMM
            #### NOTE: Pkg.add("Clp") if you don't have it

#### NOTE: Complete the missing code parts that are requested.
####       Code lines to be completed are marked as follows:
####       NOTE: Complete the following expression. (For example)

#### Function to solve linear programming problems
####
####   minimize    cᵀx
####   subject to  Ax = b
####                x ≥ 0
####
#### using ADMM. ρ is the penalty parameter.

function ADMM_linprog(c, A, b, ρ)

    t_start = time()          # Start timer
    MAX_ITER = 1000           # Max number of iterations

    (m,n) = size(A)           # Dimensions of matrix A

    x = zeros(n)              # Store x variable values
    z = zeros(n)              # Store z variable values
    u = zeros(n)              # Store u variable values

    objval = zeros(MAX_ITER)  # Store objective values
    r_norm = zeros(MAX_ITER)  # Store primal residual norms
    s_norm = zeros(MAX_ITER)  # Store dual residual norms

    #### Main ADMM iteration loop
    for k = 1:MAX_ITER

        #### NOTE: Complete this to compute x-step for the current iteration.
        ####       Compare with Exercise 9.2 to solve the KKT system but here
        ####       in this linear case instead of quadratic one. You don't need
        ####       to store the value of the Lagrangian dual variable v.
        KKT_matrix = [(eye(n)*ρ) A'; A zeros(m, m)]
        KKT_vector = -[c - ρ * (z - u); -b]
        xv = KKT_matrix \ KKT_vector
        x  = xv[1:n]

        #### NOTE: Complete this to compute z-step for the current iteration.
        ####       Compare with Exercise 9.2. Store the value of z to the
        ####       array z. Before updating z, we store the value of current z in
        ####       a local vector since it's needed when computing dual residual.
        z_old = copy(z)

        z = x + u

        for i in 1:length(z)
                if z[i] < 0
                        z[i] = 0
                end
        end

        #### NOTE: Complete this to compute u-step for the current iteration.
        ####       Compare with Exercise 9.2.
        u = u + x - z

        #### NOTE: Compute the objective value, the primal residual, and the
        ####       dual residual for the current iteration k. Put the formulas
        ####       for the residuals inside the norm function. Compare with
        ####       exercise 9.2 on how to compute the residuals.
        objval[k]  = (c' * x)[1]                     # Objective value
        r_norm[k]  = norm(x - z)             # Primal residual norm
        s_norm[k]  = ρ * norm(z - z_old)             # Dual residual norm

        #### Print progress every 5 iterations
        if k % 5 == 0
            #println("iteration: ", k ,"/\tprimal residual: ", round(r_norm[k],4),
                    #"\tdual residual: ", round(s_norm[k],4))
        end

        #### Check stopping condition.
        #### NOTE: You can try to change the tolerance 0.1
        if r_norm[k] + s_norm[k] < 0.1
            #### Cut off unnecessary values
            objval = objval[1:k]
            r_norm = r_norm[1:k]
            s_norm = s_norm[1:k]
            break
        end
    end

    #### Duration of the ADMM LP method
    t_end = time() - t_start

    return x, objval, r_norm, s_norm, t_end
end

#### LP solver Clp to compare the ADMM solution
function linprog(c, A, b)

    t_start = time()
    m,n = size(A)

    lp = Model(solver = ClpSolver())

    @variable(lp, x[1:n] >= 0)
    @constraint(lp, A*x .== b)
    @objective(lp, Min, sum(c[i]*x[i] for i = 1:n))

    solve(lp)
    t_end = time() - t_start;

    return getvalue(x), getobjectivevalue(lp), t_end
end

#### Compare ADMM with Clp on random LP problem instances
#### NOTE: You can change argument in srand(0) to generate different
####       random problem instances (e.g., srand(1))
srand(100)

n = 500                # Dimension of x
m = 400                # Number of equality constraints

c  = rand(n,1) + 0.5   # Create nonnegative price vector with mean 1
x0 = abs.(randn(n,1))  # Create random solution vector

A = abs.(randn(m,n))   # Create random, nonnegative matrix A
b = A*x0               # Create vector b

ρ = 10             #### Value of the penalty parameter to use
                       #### NOTE: You can try different values of ρ

#### Call both solvers
(x1, obj_val1, r_norm, s_norm, sol_time1) = ADMM_linprog(c, A, b, ρ)
(x2, obj_val2, sol_time2) = linprog(c, A, b)

#### Print solution times and objective values
println("\nADMM solution time: ", sol_time1[end])
println(" Clp solution time: ", sol_time2)

println("\nADMM objective value: ", obj_val1[end])
println(" Clp objective value: ", obj_val2)


j = 50
iterations = 20
rhos = linspace(0.2, 2, j)
results_ADMM = zeros(5, j)
results_clp = zeros(2, j)

for i = 1:j
        println(i)
        n = 500                # Dimension of x
        m = 400                # Number of equality constraints

        c  = rand(n,1) + 0.5   # Create nonnegative price vector with mean 1
        x0 = abs.(randn(n,1))  # Create random solution vector

        A = abs.(randn(m,n))   # Create random, nonnegative matrix A
        b = A*x0               # Create vector b

        ρ = rhos[i]             #### Value of the penalty parameter to use
                               #### NOTE: You can try different values of ρ
        result_ADMM = zeros(5)
        result_clp = 0

        for it = 1:iterations
                n = 500                # Dimension of x
                m = 400                # Number of equality constraints

                c  = rand(n,1) + 0.5   # Create nonnegative price vector with mean 1
                x0 = abs.(randn(n,1))  # Create random solution vector

                A = abs.(randn(m,n))   # Create random, nonnegative matrix A
                b = A*x0               # Create vector
                #### Call both solvers
                (x1, obj_val1, r_norm, s_norm, sol_time1) = ADMM_linprog(c, A, b, ρ)
                (x2, obj_val2, sol_time2) = linprog(c, A, b)

                result_ADMM += [length(obj_val1), (obj_val1[end]-obj_val2)^2, r_norm[end], s_norm[end], sol_time1]
                result_clp += sol_time2
        end


        results_ADMM[:,i] = result_ADMM ./ iterations
        results_clp[:,i] = result_clp / iterations
end

fig = plot(xlabel = "Penalty term", ylabel = "Time", size = (1200, 800))
plot!(rhos, results_ADMM[5,:], label = "Time")
#savefig("rho_time.pdf")

fig = plot(xlabel = "Penalty term", ylabel = "Residual", size = (1200, 800))
plot!(rhos, results_ADMM[4,:], label = "Dual residual")
plot!(rhos, results_ADMM[3,:], label = "Primal residual")
#savefig("residuals.pdf")

fig = plot(xlabel = "Penalty term", ylabel = "Iterations", size = (1200, 800))
plot!(rhos, results_ADMM[1,:], label = "Average iterations")
#savefig("iterations.pdf")

fig = plot(xlabel = "Penalty term", ylabel = "Accuracy", size = (1200, 800))
plot!(rhos, results_ADMM[2,:], label = "Average distance from optimal value")
#savefig("accuracy.pdf")

plot(rhos, results_ADMM[5,:] ./ results_clp[2,:])

plot(results_clp[1,:])
