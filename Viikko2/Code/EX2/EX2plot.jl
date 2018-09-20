## Plot the function in EX2
using PyPlot

n = 1000
x = linspace(0.5,2,n)
y = linspace(0.5,2,n)

f(x,y) = 1/(x+y)

# Compute each point in the grid (x,y) separately
z = [f(x[i],y[j]) for j = 1:n, i = 1:n]

# Plot the funcion
figure(figsize = (12,7))
surf(x,y,z)
# Optimal solution
plot3D([1],[1],[0.5],"r.",markersize = 20)
xlabel("x"); ylabel("y"); zlabel("z")
tight_layout()
