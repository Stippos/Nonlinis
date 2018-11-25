using PyPlot
using JuMP
using Ipopt

# Read the daily prices data from a .csv file.
data   = readcsv("prices.csv")
# First row has the names of the stocks
names  = data[1, 1:end-2]
# Last two columns has data and US$ rate
prices = data[2:end, 1:end-2]

#form statistics for the stocks based on the data
returns = diff(prices) ./ prices[1:end-1,:]
T, n = size(returns)      # Number of days and stocks
μ = vec(mean(returns, 1)) # Expected returns
Σ = cov(returns)          # Covariance matrix
ρ = cor(returns)          # Correlation matrix
σ = sqrt.(diag(Σ))        # Standard deviation

ix = sortperm(μ)          # ix = index ordering

#arrays for storing the results
all_x = []
all_obj = []

#loop over all the different lambdas
for lambda in 0:0.1:1
    #debugging print
    print("")
    print("LAMBDA")
    print(lambda)
    print("")

    #form the model
    m = Model(solver = IpoptSolver()) # Model

    #the wights for the portfolio
    @variable(m, x[1:n] >= 0)         # Stock positions

    #objective function as defined in the homework
    @expression(m, objective, lambda * dot(μ, x) - (1 - lambda) * dot(x, Σ*x))

    #maximize the objective
    @objective(m, Max, objective)

    #only constraint we have is that the weights must add up to one
    @constraint(m, sum(x) == 1)

    # Solve the problem and get solution
    solve(m)
    x = getvalue(x)        # Stock positions
    ret = dot(μ,x)         # Return
    std = sqrt(dot(x,Σ*x)) # Risk: std. deviation of returns
    obj = getobjectivevalue(m)

    #save the results to the vectors
    push!(all_x, x)
    push!(all_obj, obj)
end

#solve returns from the weights
returns = [dot(μ, x) for x in all_x]

#solve volatilities from the weights
vars = [dot(x, Σ*x) for x in all_x]

#plot the returns as a function of lambda
figure(figsize = (10,6))
plot(0:0.1:1, returns)
xlabel("λ")
ylabel("μ")
savefig("returns.png")

#plot the volatilities as a function of lambda
figure(figsize = (10,6))
plot(0:0.1:1, vars)
xlabel("λ")
ylabel("σ")
savefig("volatilities.png")

#plot the volatities as a function of return
figure(figsize = (10,6))
plot(returns, vars)
xlabel("μ")
ylabel("σ")
savefig("volvsret.png")
