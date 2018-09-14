#Pkg.add("JuMP")
#Pkg.add("Plots")
using Plots

## Tehtävä 1
f(x) = 2x^4 - 5x^3 - x^2
plot(f, -1, 3)

f(x) = 2x^4 - 5x^3 - x^2
plot(f, -1/7, 0)

## Tehtävä 2
Pkg.add("PyPlot")
ENV["PYTHON"]=""
Pkg.add("PyCall")
using PyCall
using PyPlot

g(x,y) = (y - x^2)^2 - x^2

x = -2:0.01:2
y = x'
z = g.(x, y)

surf(x, y, z)
