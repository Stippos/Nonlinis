#Pkg.add("Ipopt")
#Pkg.add("JuMP")

using Ipopt
using JuMP

type Arc
    i::String
    j::String
end

# type for node properties
type Node
    q1::Float64  # Propery value (sulfur %)
    q2::Float64
    ubq1::Float64
    lbq1::Float64
    ubq2::Float64
    lbq2::Float64
    flb::Float64
    fub::Float64  # Upper flow bound
end

V = ["s1", "s2", "s3", "s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","s16","p1", "p2", "t1", "t2"]
S = ["s1", "s2", "s3", "s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","s16"]
P = ["p1", "p2"]
T = ["t1", "t2"]
inode = ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10","s11","s12","s13","s14","s15","s16","p1", "p1", "p2", "p2"]
jnode = ["p1", "p1", "p1", "p1", "p1", "p1", "p1", "p1", "p2", "p2", "p2",  "p2", "p2", "p2", "p2", "p2","t1", "t2","t1", "t2"]
fubound = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,15,15,11,17]
flbound = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0 ,0 ,10,11]
lq1 = [5.84, 5.4, 0.24, 2.01, 5.85, 5.38, 0.26, 2.04, 0.64, 0.57, 0.02, 0.14, 0.93, 0.85, 0.03, 0.26, 0, 0, 0, 0]
uq1 = [5.84, 5.4, 0.24, 2.01, 5.85, 5.38, 0.26, 2.04, 0.64, 0.57, 0.02, 0.14, 0.93, 0.85, 0.03, 0.26, 100, 100, 1.5, 3.5]
lq2 = [43.7, 36.8, 12.8, 15.4, 47.3, 39.2, 13.1, 15.9, 39.9, 38.2, 13.5, 16.3, 38.1, 34.1, 13.2, 15.5, 0, 0, 30, 32]
uq2 = [43.7, 36.8, 12.8, 15.4, 47.3, 39.2, 13.1, 15.9, 39.9, 38.2, 13.5, 16.3, 38.1, 34.1, 13.2, 15.5, 100, 100, 34, 40]

A = Arc.(inode, jnode)
N = Node.(uq1, uq2, uq1, lq1, uq2, lq2, flbound, fubound)

N = Dict(V[i] => N[i] for i = 1:size(V,1))

uppper_q_bound_11 = 1.5
uppper_q_bound_12 = 3.5
uppper_q_bound_21 = 34
uppper_q_bound_22 = 40


N["p1"].q1 = 50
N["p1"].q2 = 50
N["p2"].q1 = 50
N["p2"].q2 = 50

m = Model(solver = IpoptSolver())

@variable(m, x[A] >= 0)       # Arc flows
@variable(m, q1[i in V] >= 0)
@variable(m, q2[i in V] >= 0)


## Cost and revenue
#@expression(m, revenue, sum((100 * (2 - q1[i] / N[i].ubq1) + 150 * (2 - N[i].q2 / N[i].ubq2)) * x[a] for i in T, a in A if a.j == i))

@expression(m, revenue, 100*(2 - (q1["t1"] / N["t1"].ubq1)) * (x[A[17]] + x[A[18]]) + 150*(2 - (q1["t1"] / N["t2"].ubq1)) * (x[A[19]] + x[A[20]]))

@objective(m, Max, revenue)

@constraint(m, [p in P], sum(x[a] for a in A if a.j == p) == sum(x[a] for a in A if a.i == p))
@constraint(m, [t in T], N[t].flb <= sum(x[a] for a in A if a.j == t) <= N[t].fub)

@constraint(m, [p in P], N[p].flb <= sum(x[a] for a in A if a.j == p) <= N[p].fub)

@constraint(m, [a in A], N[a.i].flb <= x[a] <= N[a.i].fub)

@constraint(m, [s in S], N[s].flb <= sum(x[a] for a in A if a.i == s) <= N[s].fub)

@constraint(m, [p in P], sum(q1[a.i]*x[a] for a in A if a.j == p) == q1[p]*sum(x[a] for a in A if a.i == p))
@constraint(m, [p in P], sum(q2[a.i]*x[a] for a in A if a.j == p) == q2[p]*sum(x[a] for a in A if a.i == p))

# Sulfur balances at T nodes
@constraint(m, [t in T], sum(q1[a.i]*x[a] for a in A if a.j == t) == q1[t]*sum(x[a] for a in A if a.j == t))
@constraint(m, [t in T], sum(q2[a.i]*x[a] for a in A if a.j == t) == q2[t]*sum(x[a] for a in A if a.j == t))

# Sulfur upper bounds at T nodes
@constraint(m, [t in T], N[t].lbq1 <= q1[t] <= N[t].ubq1)
@constraint(m, [t in T], N[t].lbq2 <= q2[t] <= N[t].ubq2)

@constraint(m, [p in P], N[p].lbq1 <= q1[p] <= N[p].ubq1)
@constraint(m, [p in P], N[p].lbq2 <= q2[p] <= N[p].ubq2)

@constraint(m, [s in S], N[s].lbq1 <= q1[s] <= N[s].ubq1)
@constraint(m, [s in S], N[s].lbq2 <= q2[s] <= N[s].ubq2)

println(m)

## Solve model
status = solve(m)
println(status)

## Get objecive + solution
obj = getobjectivevalue(m)
x   = getvalue(x)
q1   = getvalue(q1)
q2   = getvalue(q2)

println("\n\nSolution cost: ", round(obj,4))
println("\nFlows: \n")
## Arcs
for a in A
    println("x($(a.i),$(a.j)) = ", round(x[a], 4))
end

println("\nProperty values:\n")
## Property values
for p in P
    println("q1[$p] = ", round(q1[p], 4))
    println("q2[$p] = ", round(q2[p], 4))
end
for t in T
    println("q1[$t] = ", round(q1[t], 4))
    println("q2[$t] = ", round(q2[t], 4))
end
println()
