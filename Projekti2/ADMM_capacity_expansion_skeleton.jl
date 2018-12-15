addprocs()                    ###### Use max number of parallel threads
                              # NOTE: The macro @everywhere makes the loaded
                              #       packages available to all threads
@everywhere using Suppressor  ###### Suppress annoying output from each thread
@everywhere using JuMP        # NOTE: type Pkg.add("Suppressor") to install this
@everywhere using Gurobi      # Gurobi is needed for dual simplex with QP support
#### NOTE: Download the latest version of Gurobi (8.1) from
####       http://www.gurobi.com/downloads/download-center. Create an account
####       and after installing Gurobi request the academic licence which should
####       be available for Aalto students. Use Pkg.add("Gurobi") to install the
####       Gurobi solver interface for Julia.

################################################################################
## NOTE:  Use this problem instance to complete the missing parts of the code ##
##        to get your ADMM implementation working.                            ##
##################### PROBLEM DATA FOR SMALLER INSTANCE ########################
#srand(1)                          # Control random number generation
#nI = 20                           # Number of suppliers
#nJ = 30                           # Number of clients
#nS = 500                          # Number of scenarios
#I = 1:nI                          # Supplier index set
#J = 1:nJ                          # Client index set
#S = 1:nS                          # Scenario index set
#### Generate random data for the problem
#c = rand(5:20, nI)                # Unit capacity costs of suppliers
#d = rand(nJ,nS).*rand(10:50, nJ)  # Client demands in all scenarios
#q = rand(5000:10000, nJ)          # Unit costs of unfulfilled demand
#p = ones(nS).*1/nS                # Scenario probabilities
#f = rand(3:45, (nI,nJ))           # Unit costs to fulfil demands
#b = rand(20:100, nI)              # Max supplier capacities
#B = 3000                          # Max budget (cost) for capacity acquisition
################################################################################


################################################################################
## NOTE:  Use this problem instance and comment over the smaller problem data ##
##        when your algorithm works correctly.                                ##
##################### PROBLEM DATA FOR LARGER INSTANCE #########################
 srand(100)                          # Control random number generation
 nI = 20                           # Number of suppliers
 nJ = 100                          # Number of clients
 nS = 3000                         # Number of scenarios
 I = 1:nI                          # Supplier index set
 J = 1:nJ                          # Client index set
 S = 1:nS                          # Scenario index set
# #### Generate random data for the problem
 c = rand(20:100, nI)              # Unit capacity costs of suppliers
 d = rand(nJ,nS).*rand(10:50, nJ)  # Client demands in all scenarios
 q = rand(5000:10000, nJ)          # Unit costs of unfulfilled demand
 p = ones(nS).*1/nS                # Scenario probabilities
 f = rand(3:45, (nI,nJ))           # Unit costs to fulfil demands
 b = rand(200:1000, nI)            # Max supplier capacities
 B = 30000                         # Max budget (cost) for capacity acquisition
################################################################################

### We first solve the problem formulation directly without ADMM
# m = Model(solver = GurobiSolver())            # We use Gurobi solver
#
# @variable(m, x[i in I] >= 0)                  # Reserved capacity variables
# @variable(m, y[i in I, j in J, s in S] >= 0)  # Demand fulfilment variables
# @variable(m, u[j in J, s in S] >= 0)          # Unfulfilled demand variables
#
# # Objective: Minimize the total phase 1 + phase 2 costs over all scenarios
# @objective(m, Min,
#       sum(c[i]*x[i] for i in I) +
#       sum(p[s]*(f[i,j]*y[i,j,s] + q[j]*u[j,s]) for i in I, j in J, s in S))
#
# @constraint(m, [i in I], x[i] <= b[i])          # Max capacity constraint
# @constraint(m, sum(c[i]*x[i] for i in I) <= B)  # Max capacity budget (cost) constraint
# # Capacity reserve limit constraint for each supplier i in each scenario s
# @constraint(m, [i in I, s in S], sum(y[i,j,s] for j in J) <= x[i])
# # Demand balance constraint (u[j,s] is unfulfilled demand of client j in scenario s)
# @constraint(m, DemBal[j in J, s in S], sum(y[i,j,s] for i in I) + u[j,s] == d[j,s])
#
# @time solve(m)                                  # Solve the problem
# println(getvalue(x))                            # Get optimal x values (reserved capacities)
# println(getvalue(sum(c[i]*x[i] for i in I)))    # Optimal cost of reserved capacities

#rhos = [50, 100, 200, 500, 1000, 2000]
#times = [83.577, 290.998, 406.233, 780.992, 1409.77, 2698.96]
#rhos = [20, 30, 40, 50, 60, 70, 80, 90, 100, 150]
#times = [27.89, 113.739, 106.949, 104.305, 106.009, 143.728, 143.616, 142.003, 178.744, 212.476]
#seed = 10
#rhos = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#times = [28.321, 79.907, 72.589, 71.391, 106.572, 107.748, 106.522, 105.443, 144.861, 144.038, 142.769, 176.981]
#seed = 100
#rhos = [0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000]
#times = [27.602, 19.914, 20.336, 19.747, 28.945, 38.42, 47.032, 56.582, 64.782, 73.567, 99.815, 134.832, 190.175, 571.949, 304.029, 574.272]
#iterations = [2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 11, 15, 21, 28, 34, 67]
#seed = 1000
#rhos = [0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000]
#times = [52.123, 43.666, 43.573, 43.22, 52.292, 60.93, 60.396, 69.917, 80.162, 78.885, 113.55, 139.329, 201.517, 261.66, 319.523, 606.515]
#iterations = [5,5,5,5,6,7,7,8,9,9,13,16,23, 30, 37, 72]
#seed = 10000
#rhos = [0.1, 0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000, 2000]
#times = [683.599, 195.362, 196.705, 197.162, 193.141, 189.797, 190.503, 190.404, 129.177, 164.21, 190.451, 124.773, 184.352, 134.044, 194.575, 141.633, 185.023, 293.135]
#iterations = [70, 22, 22, 21, 22, 22, 22, 22, 15, 19, 22, 21, 14, 21, 14, 17, 16, 20, 32]
#seed = 100000
#rhos = [0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000]
#times = [26.638, 18.371, 17.736, 26.509, 35.934, 54.236, 68.088, 81.476, 99.267, 122.358, 162.1089, 196.334, 294.842, 382.469, 483.233, 911.7389]
#iterations = [2, 2, 2, 3, 4, 6, 7, 9, 11, 12, 17, 22, 33, 43, 54, 106]
#seed = 1000000
#rhos = [0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000]
#times = [27.404, 19.194, 18.646, 18.64, 18.568, 28.132, 37.354, 37.761, 47.092, 56.69, 83.541, 93.486, 131.828, 158.419, 193.333, 372.3619]
#iterations = [2, 2, 2, 2, 2, 3, 4, 4, 5, 6, 8, 10, 14, 17, 21, 41]
####       are requested. Code lines to be completed are marked as follows:
####       NOTE: Complete the following expression. (For example)
rhos = [0.5, 1, 5, 10, 20, 40, 50, 70, 90, 100, 150, 200, 300, 400, 500, 1000]
times = zeros(length(rhos))
for r = 1:length(rhos)
    #### Problem parameters and memory allocation. NOTE: Compare with Exercise 9.3.
    tini = time()                      # Start timer
    ρ    = rhos[r]                   # Penalty parameter. ρ ∈ [1000, 3000] should work well.
    x_s  = SharedArray(zeros(nI, nS))  # Store each x vector in (x,y)-steps of each scenario
    v_s  = SharedArray(zeros(nI, nS))  # Store each v vector in v-steps of each scenario
    z    = SharedArray(zeros(nI))      # Store z vectors at each ADNM iteration
    #### NOTE: The array type SharedArray makes the values of x_s v_s and z available
    ####       for all active threads so that they can modify and use them in parallel
    println()
    println("Rho = ", ρ)
    println()
    #### Main loop with 200 ADMM iterations
    for k = 1:200
        #### Loop for solving each (x,y) step separately for each scenario s in S
        #### NOTE: Compare with Exercise 9.3.
        #### NOTE: The macro @parallel allocates scenario subproblems to available
        ####       threads. The macro @sync ensures that all scenario subproblems
        ####       are solved before proceeding with the code after the for loop
        @sync @parallel for s in S
            #### Model to solve (x,y)-step of the current scenario subproblem s
            scen_m = Model(solver = GurobiSolver(OutputFlag = 0, Threads = 1))
            @variable(scen_m, x[i in I] >= 0)          # Reserved capacity variables
            @variable(scen_m, y[i in I, j in J] >= 0)  # Demand fulfilment variables
            @variable(scen_m, u[j in J] >= 0)          # Unfulfilled demand variables

            @constraint(scen_m, [i in I], x[i] <= b[i])          # Max capacity constraint
            @constraint(scen_m, sum(c[i]*x[i] for i in I) <= B)  # Max capacity budget (cost)
            # Capacity reserve limit constraint for each supplier i
            @constraint(scen_m, [i in I], sum(y[i,j] for j in J) <= x[i])
            # Demand balance constraint (u[j] is unfulfilled demand of client j)
            @constraint(scen_m, [j in J], sum(y[i,j] for i in I) == d[j,s] - u[j])

            #### NOTE: Complete this objective to compute (x,y)-step for the
            ####       current scenario. Compare with Exercise 9.3.

            @expression(scen_m, first, sum((c[i] + v_s[i, s]) * x[i] for i in I))
            @expression(scen_m, second, sum(f[i, j] * y[i, j] + q[j] * u[j] for i in I, j in J))
            @expression(scen_m, third, ρ/2 * sum((x[i] - z[i])^2 for i in I))

            @objective(scen_m, Min, first + second + third)

            #### Solve the (x,y) step for the current scenario
            @suppress solve(scen_m)

            #### Store the value of x for the current scenario
            x_s[:,s] = getvalue(x[:])
        end

        #### Compute primal and dual residuals. We use array to exploit parallelism
        tol = SharedArray(zeros(nS))
        #### NOTE: Complete the computation of the residual for each s.
        ####       Compare with Exercise 9.3.
        @sync @parallel for s in S
             tol[s] = p[s] * norm(x_s[:,s] - z)^2
        end
        #### Total residual = sum of subproblem residuals
        tol = sum(tol[s] for s in S)

        #### Print current progress
        println("iteration: ", k ,"/ residual: ", tol)
        #### Stopping condition: if primal + dual residual is small enough
        if tol < 1e-2
            break
        end

        #### Compute z-step for this iteration
        #### NOTE: Complete the z-step. Compare with Exercise 9.3
        z = x_s * p

        #### Update v-step separately for each scenario
        #### NOTE: Complete the v-steps for each scenario s.
        ####       Compare with exercise 9.3
        @sync @parallel for s in S
            v_s[:,s] = v_s[:,s] + ρ*(x_s[:,s] - z)
        end
    end

     #### Stop timer
     tend = time() - tini
     times[r] = tend
     #### Solution time of ADMM:
     println("\nADMM Solution time: ", tend, "\n")

     #### Kill threads

     ####  Print optimal x values (reserved capacities)
     for i = 1:nI
         println("x[$i] = ", round(x_s[i,1],2))
     end
     #### Print optimal cost of reserved capacities
     println("\nCost of reserved capacity: ", dot(c,x_s[:,1]), "\n")

end

@suppress rmprocs(procs())


print(times)
