using PyPlot       # For plotting
using ForwardDiff  # For computing gradients using automatic differentiation

# Function to minimize
f(x) = (x[1] - 2)^4 + (x[1] - 2x[2])^2

eps = 1e-04        # Tolerance for stopping criterion
N   = 100          # Number of iterations
xn  = zeros(N+1,2) # Store values xn of Newton
xn[1,:] = [0,3]    # Initial value

# Gradient and Hessian
∇(f, x)  = ForwardDiff.gradient(f, x)
∇²(f, x) = ForwardDiff.hessian(f, x)

# Newtown step
for i = 1:N
    xn[i+1,:] = xn[i,:] - inv(∇²(f, xn[i,:]))*∇(f, xn[i,:])
    # Check stopping condition
    if sum(abs.(∇(f, xn[i+1,:])) .> eps) == 0
        # Cut off unnecessary values
        xn = xn[1:i+1,:]
        break
    end
end

# Values of Newton iterations
fn = [f(xn[i,:]) for i = 1:size(xn,1)]

# Plot results
n = 1000
x1 = linspace(-1,4,n)
x2 = linspace(-1,4,n)
# Original function values
z = [f([x1[i],x2[j]]) for j = 1:n, i = 1:n]

figure(figsize = (12,8))
# Plot contours of z = f(x) and set values for which contours to plot
cs = contour(x1,x2,z, levels=[0.1;0.6;2;collect(5:10:100)])
clabel(cs, fontsize = 9)
plot(xn[:,1], xn[:,2], ".-")  # Plot Newton iterations in same figure
xlabel("x₁")
ylabel("x₂")
legend(["Newton"])
title("Newton's method")
tight_layout()
