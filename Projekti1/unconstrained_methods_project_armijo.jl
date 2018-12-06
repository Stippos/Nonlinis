using Plots, ForwardDiff, LaTeXStrings
pyplot()

include("line_searches.jl");

# function to be optimised
#f(x) = 0.26 * (x[1]^2 + x[2]^2) - 0.48 * x[1] * x[2] #uncommet for A
f(x) = exp(x[1] + 3*x[2] - 0.1) + exp(x[1] - 3*x[2] - 0.1) + exp(-x[1] - 0.1)
#f(x) = (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2 #uncomment for C


# Calculating first- and second-order derivatives
∇(f,x)  = ForwardDiff.gradient(f, x);
H(f,x) = ForwardDiff.hessian(f, x);

# Plotting the contours of the function to be optimised
n = 1000
x = linspace(-10,25,n);
y = linspace(-10,10,n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

contour(x,y,z,
        levels = [0, 0.1, 0.3, 0.5 , 1, 2, 3, 4, 6],
        xaxis = (L"$x_1$", (-1,8)),
        yaxis = (L"$x_2$", (-1,6)),
        clims = (0,20),
        clabels = true,
        aspect_ratio = :equal)

contour(x,y,z,
        levels = [2.6, 3, 4, 5 , 6, 7, 8, 10, 12],
        xaxis = (L"$x_1$", (-2,3)),
        yaxis = (L"$x_2$", (-1.5,2)),
        clims = (0,20),
        clabels = true,
        aspect_ratio = :equal)

contour(x,y,z,
        levels = [1, 5, 20, 40,70 , 100, 120,150, 170, 200],
        xaxis = (L"$x_1$", (-5,4)),
        yaxis = (L"$x_2$", (-5,4)),
        clims = (0,200),
        clabels = true,
        aspect_ratio = :equal)


gui()


#Common inititalisation parameters
M = 10001 # 10000 iter + start
ϵ = 1e-6 # tolerance
#xstart = [7;3] #uncomment for A
xstart = [1;1.5] #uncomment for B
#xstart = [-2;1] #uncomment for C

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
        λbar = newton(f,1) # λ₀ = 1 helps preventing divergence
        #λbar = golden_ratio(f,0,10)
    end
    return λbar
end


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
    λbar = line_search(ls, ARMIJO)

    xg[:,k+1] = xg[:,k] + λbar*d
    # uncomment to see progress of the algorithm
    #println("iter=", k, " λ=",λbar, " xᵏ=",xg[:,k+1])
end

tend = time() - tini; # Stop stopwatch
println("Gradient descent converged.")
println(" Total steps: ", size(xg,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xg[:,end], "/ Opt. value: ", f(xg[:,end]), "\n")

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
    λbar = line_search(ls, ARMIJO)

    xn[:,k+1] = xn[:,k] + λbar*d
    # Uncomment to see progress of the algorithm
    #println("iter.=", k, " λ=",λbar, "xᵏ=",xg[:,k+1])
end
tend = time() - tini; # Stop stopwatch
println("Newton's method converged.")
println(" Total steps: ", size(xn,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xn[:,end], "/ Opt. value: ", f(xn[:,end]), "\n")

#Plotting progress of Newton's method
plot!(xn[1,:], xn[2,:],label = "Newton", marker=:circle)


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
    λbar = line_search(ls, ARMIJO)

    # Update function difference and do a step
    p = λbar*d
    xb[:,k+1] = xb[:,k] + p

    # Check stop criterion using the norm of ∇f(x)
    if norm(∇(f, xb[:,k+1])) < ϵ
        xb = xb[:,1:k+1]
        break
    end

    # Update Gradient difference
    q   = ∇(f, xb[:,k+1]) - ∇(f, xb[:,k])
    # Update Hessian approximation
    B = B + (q*q')/(q'*p) - (B*p*p'*B)/(p'*B*p)
end
tend = time() - tini; # Stop stopwatch
println("Quasi-Newton (BFGS) converged.")
println(" Total steps: ", size(xb,2)-1)
println(" Total time (s): ", tend)
println(" Sol. found: ", xb[:,end], "/ Opt. value: ", f(xb[:,end]), "\n")

plot!( xb[1,:], xb[2,:], label = "BFGS", marker=:circle)


#savefig("C_paths_armijo.pdf")
# Comparing convergence between algorithms
dist_xg = sqrt.(sum(( xg .- [0, 0]).^2,1)');
dist_xn = sqrt.(sum(( xn .- [0, 0]).^2,1)');
dist_xb = sqrt.(sum(( xb .- [0, 0]).^2,1)');

dist_xg = sqrt.(sum(( xg .- [-0.346574, 0]).^2,1)');
dist_xn = sqrt.(sum(( xn .- [-0.346574, 0]).^2,1)');
dist_xb = sqrt.(sum(( xb .- [-0.346574, 0]).^2,1)');

dist_xg = sqrt.(sum(( xg .- [-2.80512, 3.13131]).^2,1)');
dist_xn = sqrt.(sum(( xn .- [-3.77931, -3.28319]).^2,1)');
dist_xb = sqrt.(sum(( xb .- [-3.77931, -3.28319]).^2,1)');

plot(dist_xg, yscale=:log10, label = "Gradient")
plot!(dist_xn, yscale=:log10, label = "Newton")
plot!(dist_xb, yscale=:log10, label = "BFGS",
    xaxis = ("iterations", (1,50)),
    yaxis = (L"$||x_k - \overline{x}||$", ( ϵ, 10)))

#savefig("C_convergence_armijo.pdf")
