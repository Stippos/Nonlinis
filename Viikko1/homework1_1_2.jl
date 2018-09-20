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

V = ["s1", "s2", "s3", "s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","s16", "t1", "t2"]
S = ["s1", "s2", "s3", "s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","s16"]
T = ["t1", "t2"]
inode = ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10","s11","s12","s13","s14","s15","s16", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10","s11","s12","s13","s14","s15","s16"]
jnode = ["t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2", "t2"]
fubound = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,11,17]
flbound = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0 ,10,11]
lq1 = [5.84, 5.4, 0.24, 2.01, 5.85, 5.38, 0.26, 2.04, 0.64, 0.57, 0.02, 0.14, 0.93, 0.85, 0.03, 0.26, 0, 0]
uq1 = [5.84, 5.4, 0.24, 2.01, 5.85, 5.38, 0.26, 2.04, 0.64, 0.57, 0.02, 0.14, 0.93, 0.85, 0.03, 0.26, 1.5, 3.5]
lq2 = [43.7, 36.8, 12.8, 15.4, 47.3, 39.2, 13.1, 15.9, 39.9, 38.2, 13.5, 16.3, 38.1, 34.1, 13.2, 15.5, 30, 32]
uq2 = [43.7, 36.8, 12.8, 15.4, 47.3, 39.2, 13.1, 15.9, 39.9, 38.2, 13.5, 16.3, 38.1, 34.1, 13.2, 15.5, 34, 40]

A = Arc.(inode, jnode)
N = Node.(uq1, uq2, uq1, lq1, uq2, lq2, flbound, fubound)

N = Dict(V[i] => N[i] for i = 1:size(V,1))

##

m = Model(solver = IpoptSolver())

@variable(m, x[A] >= 0)       # Arc flows
@variable(m, q1[i in V] >= 0)
@variable(m, q2[i in V] >= 0)

#setvalue(q1["p1"], 1)
#setvalue(q1["p2"], 33)
#setvalue(q2["p1"], 1)
#setvalue(q2["p2"], 33)

#setvalue(q1["p1"], 0)
#setvalue(q1["p2"], 0)
#setvalue(q2["p1"], 0)
#setvalue(q2["p2"], 0)



#@expression(m, revenue, sum((100 * (2 - q1[i] / N[i].ubq1) + 150 * (2 - N[i].q2 / N[i].ubq2)) * x[a] for i in T, a in A if a.j == i))

@expression(m, revenue, 100*(2 - (q1["t1"] / N["t1"].ubq1)) * (sum(x[a] for a in A if a.j == "t1")) + ((sum(x[a] for a in A if a.j == "t2")) * 150*(2 - (q1["t2"] / N["t2"].ubq1))))

@objective(m, Max, revenue)

@constraint(m, [t in T], N[t].flb <= sum(x[a] for a in A if a.j == t) <= N[t].fub)

@constraint(m, [a in A], N[a.i].flb <= x[a] <= N[a.i].fub)

@constraint(m, [s in S], N[s].flb <= sum(x[a] for a in A if a.i == s) <= N[s].fub)

# Sulfur balances at T nodes
@constraint(m, [t in T], sum(q1[a.i]*x[a] for a in A if a.j == t) == q1[t]*sum(x[a] for a in A if a.j == t))
@constraint(m, [t in T], sum(q2[a.i]*x[a] for a in A if a.j == t) == q2[t]*sum(x[a] for a in A if a.j == t))

# Sulfur upper bounds at T nodes
@constraint(m, [t in T], N[t].lbq1 <= q1[t] <= N[t].ubq1)
@constraint(m, [t in T], N[t].lbq2 <= q2[t] <= N[t].ubq2)

@constraint(m, [s in S], N[s].lbq1 <= q1[s] <= N[s].ubq1)
@constraint(m, [s in S], N[s].lbq2 <= q2[s] <= N[s].ubq2)

#@constraint(m, [p in P], q1[p] == N[p].q1)
#@constraint(m, [p in P], q2[p] == N[p].q2)

println(m)

## Solve model
status = solve(m)
println(status)

## Get objecive + solution
obj = getobjectivevalue(m)
x_value   = getvalue(x)
q1_value   = getvalue(q1)
q2_value   = getvalue(q2)

println("\n\nSolution cost: ", round(obj,4))
println("\nFlows: \n")
## Arcs
for a in A
    println("x($(a.i),$(a.j)) = ", round(x_value[a], 4))
end

println("\nProperty values:\n")
## Property values

for t in T
    println("q1[$t] = ", round(q1_value[t], 4))
    println("q2[$t] = ", round(q2_value[t], 4))
end

println()

println("Flow to t1 = ", round(sum(x_value[a] for a in A if a.j == "t1")))
println("Flow to t2 = ", round(sum(x_value[a] for a in A if a.j == "t2")))



##
