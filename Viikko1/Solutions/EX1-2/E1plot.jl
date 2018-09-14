using PyPlot                 # Package used for plotting
pygui(true)

n = 1000                     # We use 1000 samples
x = linspace(-3,3,n)         # between [-3,3]

f(x) = 2x.^4 - 5x.^3 - x.^2  # Notice the .^ notation to handle vector inputs

figure(figsize = (10,6))     # Open figure and plot
plot(x,f(x))
xlabel("x")
ylabel("f(x)")
## Uncomment to zoom in on different parts of the plot
# axis([-0.25,0,-0.01,0.01])      # min
# axis([-0.1,0.12,-0.016,0.001])  # max
# axis([1.4,2.6,-13,-6])          # min (global)
tight_layout()
