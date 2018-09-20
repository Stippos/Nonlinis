## Exercise 2.2 + How to set inital values
## for variables before solving the problem.
using JuMP
using Ipopt

## Model and solver
m = Model(solver = IpoptSolver())
## Variables
@variable(m, x >= 0)
@variable(m, y >= 0)

## Objective (NOTE: @NLobjective instead of @objective)
@NLobjective(m, Max, 1/(x+y))

## Constraints
@constraint(m, x*y >= 1)

## Set different starting values
setvalue(x, 1000)
setvalue(y, 1000)

## Solve and print solution
solve(m)

z = getobjectivevalue(m)
x = getvalue(x)
y = getvalue(y)

println("\nCost: ", z)
println("x = $x")
println("y = $y")
