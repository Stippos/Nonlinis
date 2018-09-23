using PyPlot
using JuMP
using Ipopt

# Read the daily prices data from a .csv file.
data   = readcsv("prices.csv")
# First row has the names of the stocks
names  = data[1, 1:end-2]
# Last two columns has data and US$ rate
prices = data[2:end, 1:end-2]

returns = diff(prices) ./ prices[1:end-1,:]
T, n = size(returns)      # Number of days and stocks
μ = vec(mean(returns, 1)) # Expected returns
Σ = cov(returns)          # Covariance matrix
ρ = cor(returns)          # Correlation matrix
σ = sqrt.(diag(Σ))        # Standard deviation

ix = sortperm(μ)          # ix = index ordering

all_x = []
all_obj = []

for lambda in 0:0.1:1
    print("")
    print("LAMBDA")
    print(lambda)
    print("")
    m = Model(solver = IpoptSolver()) # Model
    @variable(m, x[1:n] >= 0)         # Stock positions

    @expression(m, objective, lambda * dot(μ, x) - (1 - lambda) * dot(x, Σ*x))

    @objective(m, Max, objective)    # Minimize variance

    @constraint(m, sum(x) == 1)
    # Solve the problem and get solution
    solve(m)
    x = getvalue(x)        # Stock positions
    ret = dot(μ,x)         # Return
    std = sqrt(dot(x,Σ*x)) # Risk: std. deviation of returns
    obj = getobjectivevalue(m)

    push!(all_x, x)
    push!(all_obj, obj)
end

returns = [dot(μ, x) for x in all_x]
vars = [dot(x, Σ*x) for x in all_x]

figure(figsize = (10,6))
plot(0:0.1:1, returns)
xlabel("λ")
ylabel("μ")
savefig("returns.png")


figure(figsize = (10,6))
plot(0:0.1:1, vars)
xlabel("λ")
ylabel("σ")
savefig("volatilities.png")

figure(figsize = (10,6))
plot(vars, returns)
xlabel("μ")
ylabel("σ")
savefig("volvsret.png")
