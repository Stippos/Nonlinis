using PyPlot        # For plotting
using ForwardDiff   # For computing derivatives with automatic differentiation

# Different Function to minimize
f(x) = x^4 - x^3 - 8x^2      # Test different initial points (Normal)
# f(x) = -x^2                #                               (Max)
# f(x) = x^5 - 6x^3 - 2x^2   # xn[1] = 1 axis([-3,4,-30,20]) (Saddle)
# f(x) = atan(x)             # xn[1] = 1 axis([-3,10,-3,3])  (Divergence)
# f(x) = (1/4)x^4 - x^2 + 2x # xn[1] = 0 axis([-3,4,-5,4])     (Cycle)

eps = 1e-06      # Tolerance for stopping criterion
N   = 100        # Number of iterations
xn  = zeros(N+1) # Store values xn of Newton
xn[1] = 1    # Initial value


# First and decond order derivatives
D(f, x)  = ForwardDiff.derivative(f, x)
D²(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)

# Newtown step
for i = 1:N
    xn[i+1] = xn[i] - D(f, xn[i])/D²(f, xn[i])
    # Check stopping condition
    if abs(D(f, xn[i+1])) < eps
        # Cut off unnecessary values
        xn = xn[1:i+1]
        break
    end
end

# Values of Newton iterations
fn = f.(xn)

# Plot results
n = 1000
x = linspace(-3,4,n)
# Original function values
z = f.(x)

figure(figsize = (12,8))
plot(x, z)                  # Plot original function
plot(xn, fn, ".-")          # Plot Newton iterations in the same figure
legend(["f(x)", "Newton"])
xlabel("x")
ylabel("f(x)")
title("Newton's method")
axis([-3,4,-30,20])         # Set axis accordingly
# axis([-3,4,-5,4])         # (Cycle)
# axis([-3,10,-3,3])        # (Divergence) (arctan)
tight_layout()
