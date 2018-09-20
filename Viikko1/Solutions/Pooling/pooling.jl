## Pooling problem with 2 start nodes, 1 pool and 2 targets
#
#  S1----\
#         P---------M
#         |
#  S2----/ \-------\
#                   M
#  S3--------------/
#
##


using JuMP          # JuMP: Modeling language and solver interface
using Ipopt         # Nonlinear programming solver
using CSV           # For reading CSV files
using DataFrames    # For arranging data to a nice format

# type for arcs (i,j)
type Arc
    i::String
    j::String
end

# type for node properties
type Node
    c::Float64  # Cost
    q::Float64  # Propery value (sulfur %)
    b::Float64  # Upper flow bound
end

## Read problem data
# Define column types
coltype = [fill(Union{String, Missing}, 6);
           fill(Union{Float64, Missing}, 3)]
# Read problem data. NOTE: Julia console must be in same folder as data.csv
# Set, e.g., Julia -> Working Directory -> Current File's Folder
data = CSV.read("data.csv", types = coltype)
println(data)                      # Print data
V = collect(skipmissing(data[1]))  # skipmissing ignores missing values
S = collect(skipmissing(data[2]))  # collect turns them into a normal array
P = collect(skipmissing(data[3]))  # NOTE: V = node labels
T = collect(skipmissing(data[4]))
inode  = collect(skipmissing(data[5]))
jnode  = collect(skipmissing(data[6]))
cost   = collect(skipmissing(data[7]))
qvalue = collect(skipmissing(data[8]))
demand = collect(skipmissing(data[9]))
# Construct arc and node types
A = Arc.(inode, jnode)
N = Node.(cost, qvalue, demand)
# Put nodes in Dictionary
N = Dict(V[i] => N[i] for i = 1:size(V,1))

# Set initial property values for nodes in P (change for different results)
N["p1"].q = 1

## Define solver
m = Model(solver = IpoptSolver())

## Variables
@variable(m, x[A] >= 0)       # Arc flows
@variable(m, q[i in V] >= 0)  # Property values

## Cost and revenue
@expression(m, cost,   sum(N[s].c*x[a] for s in S, a in A if a.i == s))
@expression(m, revenue, sum(N[t].c*x[a] for t in T, a in A if a.j == t))

## Objective
@objective(m, Max, revenue - cost)

## Constraints
# Flow balance
@constraint(m, [p in P], sum(x[a] for a in A if a.j == p) ==
                         sum(x[a] for a in A if a.i == p))
# Upper flow bound at T nodes
@constraint(m, [t in T], sum(x[a] for a in A if a.j == t) <= N[t].b)
# Sulfur balance at P nodes
@constraint(m, [p in P], sum(q[a.i]*x[a] for a in A if a.j == p) ==
                         q[p]*sum(x[a] for a in A if a.i == p))
# Sulfur balances at T nodes
@constraint(m, [t in T], sum(q[a.i]*x[a] for a in A if a.j == t) ==
                         q[t]*sum(x[a] for a in A if a.j == t))
# Sulfur upper bounds at T nodes
@constraint(m, [t in T], q[t] <= N[t].q)
# Sulfur values at S nodes
@constraint(m, [s in S], q[s] == N[s].q)
# Set initial sulfur values at P nodes
@constraint(m, [p in P], q[p] == N[p].q)

# Print model. NOTE: Print model at any point to see how it looks
println(m)

## Solve model
status = solve(m)
println(status)

## Get objecive + solution
obj = getobjectivevalue(m)
x   = getvalue(x)
q   = getvalue(q)

println("\n\nSolution cost: ", round(obj,4))
println("\nFlows: \n")
## Arcs
for a in A
    println("x($(a.i),$(a.j)) = ", round(x[a], 4))
end

println("\nProperty values:\n")
## Property values
for p in P
    println("q[$p] = ", round(q[p], 4))
end
for t in T
    println("q[$t] = ", round(q[t], 4))
end
println()
