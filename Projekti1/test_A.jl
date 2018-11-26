include("solution_methods.jl")
using ForwardDiff

start = [7, 3]

g(x) = 0.26 * (x[1]^2 + x[2]^2) - 0.48 * x[1] * x[2]
ϵ = 1e-8
iters = 10000

gd_arm_start = time();
sol_gd_arm = Gradient(g, start, iters, ARMIJO, ϵ)
gd_arm_end = time();

println("Finished gradient descent with gd_armijo in: ", sol_gd_arm[3], " steps.")
println("Time taken: ", gd_arm_end - gd_arm_start)
println("Optimal point: ", sol_gd_arm[1])
println("Optimal value: ", sol_gd_arm[2])
println("")

gd_newt_start = time();
sol_gd_newt = Gradient(g, start, iters, NEWTON, ϵ)
gd_newt_end = time();

println("Finished gradient descent with gd_newton in: ", sol_gd_newt[3], " steps.")
println("Time taken: ", gd_newt_end - gd_newt_start)
println("Optimal point: ", sol_gd_newt[1])
println("Optimal value: ", sol_gd_newt[2])
println("")

newt_newt_start = time();
sol_newt_newt = Newton(g, start, iters, NEWTON, ϵ)
newt_newt_end = time();

println("Finished Newton with newton in: ", sol_newt_newt[3], " steps.")
println("Time taken: ", newt_newt_end - newt_newt_start)
println("Optimal point: ", sol_newt_newt[1])
println("Optimal value: ", sol_newt_newt[2])
println("")

newt_arm_start = time();
sol_newt_arm = Newton(g, start, iters, ARMIJO, ϵ)
newt_arm_end = time();

println("Finished Newton with armijo in: ", sol_newt_arm[3], " steps.")
println("Time taken: ", newt_arm_end - newt_arm_start)
println("Optimal point: ", sol_newt_arm[1])
println("Optimal value: ", sol_newt_arm[2])
println("")

BFGS_newt_start = time();
sol_BFGS_newt = BFGS(g, start, iters, NEWTON, ϵ)
BFGS_newt_end = time();

println("Finished BFGS with newton in: ", sol_BFGS_newt[3], " steps.")
println("Time taken: ", BFGS_newt_end - BFGS_newt_start)
println("Optimal point: ", sol_BFGS_newt[1])
println("Optimal value: ", sol_BFGS_newt[2])
println("")

BFGS_arm_start = time();
sol_BFGS_arm = BFGS(g, start, iters, ARMIJO, ϵ)
BFGS_arm_end = time();

println("Finished BFGS with armijo in: ", sol_BFGS_arm[3], " steps.")
println("Time taken: ", BFGS_arm_end - BFGS_arm_start)
println("Optimal point: ", sol_BFGS_arm[1])
println("Optimal value: ", sol_BFGS_arm[2])
println("")
