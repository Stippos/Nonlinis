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
 srand(1)                          # Control random number generation
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

#### We first solve the problem formulation directly without ADMM
m = Model(solver = GurobiSolver())            # We use Gurobi solver

@variable(m, x[i in I] >= 0)                  # Reserved capacity variables
@variable(m, y[i in I, j in J, s in S] >= 0)  # Demand fulfilment variables
@variable(m, u[j in J, s in S] >= 0)          # Unfulfilled demand variables

# Objective: Minimize the total phase 1 + phase 2 costs over all scenarios
@objective(m, Min,
      sum(c[i]*x[i] for i in I) +
      sum(p[s]*(f[i,j]*y[i,j,s] + q[j]*u[j,s]) for i in I, j in J, s in S))

@constraint(m, [i in I], x[i] <= b[i])          # Max capacity constraint
@constraint(m, sum(c[i]*x[i] for i in I) <= B)  # Max capacity budget (cost) constraint
# Capacity reserve limit constraint for each supplier i in each scenario s
@constraint(m, [i in I, s in S], sum(y[i,j,s] for j in J) <= x[i])
# Demand balance constraint (u[j,s] is unfulfilled demand of client j in scenario s)
@constraint(m, DemBal[j in J, s in S], sum(y[i,j,s] for i in I) + u[j,s] == d[j,s])

@time solve(m)                                  # Solve the problem
println(getvalue(x))                            # Get optimal x values (reserved capacities)
println(getvalue(sum(c[i]*x[i] for i in I)))    # Optimal cost of reserved capacities


#### NOTE: ADMM approach starts from here. Complete the missing code parts that
####       are requested. Code lines to be completed are marked as follows:
####       NOTE: Complete the following expression. (For example)

#### Problem parameters and memory allocation. NOTE: Compare with Exercise 9.3.
tini = time()                      # Start timer
ρ    = 1000                        # Penalty parameter. ρ ∈ [1000, 3000] should work well.
x_s  = SharedArray(zeros(nI, nS))  # Store each x vector in (x,y)-steps of each scenario
v_s  = SharedArray(zeros(nI, nS))  # Store each v vector in v-steps of each scenario
z    = SharedArray(zeros(nI))      # Store z vectors at each ADNM iteration
#### NOTE: The array type SharedArray makes the values of x_s v_s and z available
####       for all active threads so that they can modify and use them in parallel

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
        return
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

#### Solution time of ADMM:
println("\nADMM Solution time: ", tend, "\n")

#### Kill threads
@suppress rmprocs(procs())

####  Print optimal x values (reserved capacities)
for i = 1:nI
    println("x[$i] = ", round(x_s[i,1],2))
end
#### Print optimal cost of reserved capacities
println("\nCost of reserved capacity: ", dot(c,x_s[:,1]), "\n")
