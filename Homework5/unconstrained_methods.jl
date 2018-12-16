using Plots, ForwardDiff, LaTeXStrings
pyplot()

include("line_searches.jl");

# function to be optimised
f(x) = exp(x[1] + 3*x[2] - 0.1) + exp(x[1] - 3*x[2] - 0.1) + exp(-x[1] - 0.1)

# Calculating first- and second-order derivatives
∇(f,x)  = ForwardDiff.gradient(f, x);
H(f,x) = ForwardDiff.hessian(f, x);

# Plotting the contours of the function to be optimised
n = 1000
x = linspace(-10,25,n);
y = linspace(-10,10,n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

contour(x,y,z,
        levels = [3.6, 4, 5, 6 , 8, 10, 15, 20, 30],
        xaxis = (L"$x_1$", (-5,22)),
        yaxis = (L"$x_2$", (-10,10)),
        clims = (0,20),
        clabels = true,
        aspect_ratio = :equal)

gui()

# The optimal value of x:
xopt = [-(5/6)*(-3 + 2*log(2) - 2*log(5)); 0]
fopt = f(xopt)


#Common inititalisation parameters
M = 10001 # 10000 iter + start
ϵ = 1e-6 # tolerance
xstart = [8;-3] # starting point
n = length(xstart) # dimension of the problem

# This helper function allows uschanging the line search method used
# Creating control parameter for line searches
@enum LS EXACT ARMIJO  #i.e., LS ∈ {ARMIJO, EXACT}
function line_search(f, LS)
    λbar = 0
    if LS == ARMIJO
        λbar = armijo(f)
    else # LS == EXACT
        #
        λbar = newton(f,0) # λ₀ = 1 helps preventing divergence
        #λbar = golden_ratio(f,0,10)
    end
    return λbar
end



# Method 1: coordinate descent

xc = zeros(n,M) # store step history
xc[:,1] = xstart # initial point
k = 1; # counter

println("Starting coordinate descent...")
tini = time(); # Start stopwatch
while k <= M-1
    # stop criteria on the norm of the difference between 2 sucessive points
    if   k > 1 && norm(xc[:,k] - xc[:,k-1]) < ϵ
         xc = xc[:,1:k]
         break
    end
    for j = 1:n
        d = zeros(n) # setting search direction as dⱼ = 1.
        d[j] = 1

        ls(λ) = f(xc[:,k] + λ*d) # state line search (ls) function

        λbar = golden_ratio(ls, -10, 10) # To be fully derivative independent, we use golden ratio
                                         # method. Notice it requires λ ∈ R
        xc[:,k+1] = xc[:,k] + λbar*d

        # Uncomment to see progress of the algorithm
        #println("iter=", k, " λ=",λbar, " xᵏ=",xc[:,k+1])

        k = k+1
    end
end
tend = time() - tini; # Stop stopwatch
println("Coordinate descent converged.")
println(" Total steps: ", k-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xc[:,k], "/ Opt. value: ", f(xc[:,k]), "\n")

#Plotting progress gradient descent
plot!( xc[1,:], xc[2,:],label = "Coordinate", marker=:circle)



# Method 2: gradient descent

xg = zeros(n,M) # store step history
xg[:,1] = xstart # initial point

println("Starting gradient descent...")
tini = time(); # Start stopwatch
for k = 1:M-1
    d = -∇(f,xg[:,k])
    # stop criteria on the norm of ∇f(x)
    if norm(d) < ϵ
         xg = xg[:,1:k]
         break
    end
    # state line search (ls) function
    ls(λ) = f(xg[:,k] + λ*d)

    # Line search according to selected method
    #λbar = line_search(ls, ARMIJO)
    λbar = line_search(ls, EXACT)

    xg[:,k+1] = xg[:,k] + λbar*d
    # uncomment to see progress of the algorithm
    #println("iter=", k, " λ=",λbar, " xᵏ=",xg[:,k+1])
end

tend = time() - tini; # Stop stopwatch
println("Gradient descent converged.")
println(" Total steps: ", size(xg,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xg[:,k], "/ Opt. value: ", f(xg[:,k]), "\n")

#Plotting progress gradient descent
plot!( xg[1,:], xg[2,:],label = "Gradient", marker=:circle)


### Lecture 8 methods

# Method 3: Newton's method

xn = zeros(n,M) # store step history
xn[:,1] = xstart # initial point

println("Starting Newton's method...")
tini = time(); # Start stopwatch
for k = 1:M-1
    d = -H(f,xn[:,k])\∇(f,xn[:,k])
    #stop criteria on the norm of ∇f(x)
    if norm(d) < ϵ
         xn= xn[:,1:k]
         break
    end
    #state line search (ls) function
    ls(λ) = f(xn[:,k] + λ*d)

    # Line search according to selected method
    #λbar = line_search(ls, ARMIJO)
    λbar = line_search(ls, EXACT)

    xn[:,k+1] = xn[:,k] + λbar*d
    # Uncomment to see progress of the algorithm
    #println("iter.=", k, " λ=",λbar, "xᵏ=",xg[:,k+1])
end
tend = time() - tini; # Stop stopwatch
println("Newton's method converged.")
println(" Total steps: ", size(xn,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xn[:,k], "/ Opt. value: ", f(xn[:,k]), "\n")

#Plotting progress of Newton's method
plot!(xn[1,:], xn[2,:],label = "Newton", marker=:circle)



#Method 4: Conjugate gradient method

xj = zeros(n,M)
xj[:,1] = [1, 1.5]
α = 0
d = -∇(f,xstart)
k = 1

println("Starting conjugate gradient method...")
tini = time(); # Start stopwatch
while k <= M-1
    for j = 1:length(xstart) #n=2
        #state line search (ls) function
        ls(λ) = f(xj[:,k] + λ*d)

        # Line search according to selected method
        #λbar = line_search(ls, ARMIJO)
        λbar = line_search(ls, EXACT)

        xj[:,k+1] = xj[:,k] + λbar*d
        # uncomment to see progress of the algorithm
        #println("iter=", k, " λ=",λbar, " xᵏ=",xj[:,k+1])

        #α-update following Fletcher-Reeves
        α = norm(∇(f,xj[:,k+1]))^2/norm(∇(f,xj[:,k]))^2
        d = -∇(f,xj[:,k+1]) + α*d
        k = k+1
    end
    d = -∇(f,xj[:,k])
    # stop criteria on the norm of ∇f(x)
    if norm(d) < ϵ
         xj = xj[:,1:k]
         break
    end
end
tend = time() - tini; # Stop stopwatch
println("Conjugate gradient converged.")
println(" Total steps: ", size(xj,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xj[:,k], "/ Opt. value: ", f(xj[:,k]), "\n")

#Plotting progress of conjugate gradient method
plot!( xj[1,:], xj[2,:], label = "Conjugate", marker=:circle)



#Method 5: Quasi-Newton method (BFGS)

xb = zeros(n,M)
xb[:,1] = xstart
B   = H(f,xb[:,1])
∇fᵢ = ∇(f, xb[:,1])

println("Starting quasi-Newton (BFGS) method...")
tini = time(); # Start stopwatch
for k = 1:M-1
    d = -inv(B)*∇(f, xb[:,k])

    ls(λ) = f(xb[:,k] + λ*d)

    # Line search according to selected method
    #λbar = line_search(ls, ARMIJO)
    λbar = line_search(ls, EXACT)

    # Update function difference and do a step
    p = λbar*d
    xb[:,k+1] = xb[:,k] + p

    # Check stop criterion using the norm of ∇f(x)
    if norm(∇(f, xb[:,k+1])) < ϵ
        xb = xb[:,1:k+1]
        break
    end

    # Update Gradient difference
    q  = ∇(f, xb[:,k+1]) - ∇(f, xb[:,k])
    # Update Hessian approximation
    B = B + (q*q')/(q'*p) - (B*p*p'*B)/(p'*B*p)
end
tend = time() - tini; # Stop stopwatch
println("Quasi-Newton (BFGS) converged.")
println(" Total steps: ", size(xb,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xb[:,k], "/ Opt. value: ", f(xb[:,k]), "\n")

plot!( xb[1,:], xb[2,:], label = "BFGS", marker=:circle)




# Comparing convergence between algorithms
dist_xg = sqrt.(sum(( xg .- xopt).^2,1)');
dist_xn = sqrt.(sum(( xn .- xopt).^2,1)');
dist_xj = sqrt.(sum(( xj .- xopt).^2,1)');
dist_xb = sqrt.(sum(( xb .- xopt).^2,1)');

plot(dist_xg, yscale=:log10, label = "Gradient")
plot!(dist_xn, yscale=:log10, label = "Newton")
plot!(dist_xj, yscale=:log10, label = "Conjugate")
plot!(dist_xb, yscale=:log10, label = "BFGS",
    xaxis = ("iterations", (1,15)),
    yaxis = (L"$||x_k - \overline{x}||$", ( ϵ, 10)))

# savefig("ALG_convergence.pdf")
