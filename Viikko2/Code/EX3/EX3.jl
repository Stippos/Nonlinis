using PyPlot
using JuMP
using Ipopt

# Read the daily prices data from a .csv file.
data   = readcsv("prices.csv")
# First row has the names of the stocks
names  = data[1, 1:end-2]
# Last two columns has data and US$ rate
prices = data[2:end, 1:end-2]

## Plot the data + add labels
figure(figsize = (12,6))
plot(prices); legend(names); grid("on")
xlabel("Days"); ylabel("Prices")
tight_layout()

# Returns are calculated as (p(t+1) - p(t))/p(t))*100 %
returns = diff(prices) ./ prices[1:end-1,:]
T, n = size(returns)      # Number of days and stocks
μ = vec(mean(returns, 1)) # Expected returns
Σ = cov(returns)          # Covariance matrix
ρ = cor(returns)          # Correlation matrix
σ = sqrt.(diag(Σ))        # Standard deviation

# Sort stocks by expected return:
ix = sortperm(μ)          # ix = index ordering
## Plot standard deviation (sorted order)
figure(figsize = (12,6))
subplot(211)
title("Std. deviation σ and expected return μ of 19 stocks")
x = 1:1:19; xlim(1,19)
plot(x, σ[ix], "g.-")    # Std. deviation
ylabel("σ")
tight_layout()
## Plot expected return (sorted order)
subplot(212)
x = 1:1:19; xlim(1,19)
plot(x, μ[ix], ".-")     # Expected return
plot([1,19],[0,0],"r--") # Draw line at zero
ylabel("μ")
tight_layout()

## Plot correlation matrix (sorted order)
figure(figsize = (12,6))
imshow(ρ[ix,ix], extent = [1,19,19,1]);
colorbar(); axis("image")
title("Correlation matrix of 19 stocks")
tight_layout()

## Plot each stock individually
figure(figsize = (12,6))
plot(σ[ix], μ[ix], "b.", markersize = 12)
xlabel("std deviation")
ylabel("expected return")
grid("on"); tight_layout()

## Compute the minimum risk portfolio with different average returns
μ_min = 0.0003                    # Minimum expected average return
m = Model(solver = IpoptSolver()) # Model
@variable(m, x[1:n] >= 0)         # Stock positions
@objective(m, Min, dot(x,Σ*x))    # Minimize variance
@constraint(m, dot(x,μ) >= μ_min) # Expected average return bound
@constraint(m, sum(x) == 1)
# Solve the problem and get solution
solve(m)
x = getvalue(x)        # Stock positions
ret = dot(μ,x)         # Return
std = sqrt(dot(x,Σ*x)) # Risk: std. deviation of returns

# Plot optimal asset selection
figure(figsize = (12,6))
xlim(1.5,19.5)
bar(1:19, x[ix])
title(string("Optimal assets for μ = ", round(ret,5), ", σ = ", round(std,5)))
tight_layout()

## Compute optimal tradeoff curve
N = 50
ret = zeros(N)
std = zeros(N)
μ_min_values = linspace(0,0.000869,N)

for (i,μ_min) in enumerate(μ_min_values)
    m = Model(solver = IpoptSolver(print_level = 0)) # Model
    @variable(m, x[1:n] >= 0)         # Stock positions
    @objective(m, Min, dot(x,Σ*x))    # Minimize variance
    @constraint(m, dot(x,μ) >= μ_min) # Expected average return bound
    @constraint(m, sum(x) == 1)
    # Solve the problem and get solution
    solve(m)
    x = getvalue(x)           # Stock positions
    ret[i] = dot(μ,x)         # Return
    std[i] = sqrt(dot(x,Σ*x)) # Risk: std. deviation of returns
end

## Plot tradeoff curve
figure(figsize=(12,8))
plot(std, ret, "b-")
plot(σ, μ, "k.", markersize = 12)
xlabel("std deviation")
ylabel("expected return")
grid("on"); tight_layout()

print(x)
