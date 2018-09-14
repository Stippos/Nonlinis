using PyPlot    # Package for plotting
using3D()       # We need this when plotting 3d images (e.g., surf(x,y,z))

n = 1000              # We use a 1000*1000 grid in (x,y)-plane
x = linspace(-0.5,0.5,n)  # between (x,y) ∈ ([-2,2],[-1,1])
y = linspace(-0.5,0.5,n)

f(x,y) = (y - x^2)^2 - x^2

# Here We compute each point in the grid separately
z = [f(x[i],y[j]) for j = 1:n, i = 1:n]
figure(figsize = (8,7))

# Create 2x1 subplot matrix (21) and plot the surface of f(x) in subplot 1
subplot(211, projection = "3d")
surf(x,y,z)
xlabel("x")
ylabel("y")
zlabel("z")

# Plot the contours of f(x) in (x,y)-plane in subplot 2
subplot(212)
cs = contour(x,y,z)                              #
clabel(cs, fontsize = 9, levels=[0,0.1,0.2,0.3]) # This will display contour values
xlabel("x")
ylabel("y")
tight_layout()
