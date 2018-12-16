using Plots, ForwardDiff, LaTeXStrings
pyplot()

function cjg(f, start, iterations, tolerance, newton_start)
    log = zeros(2, iterations)
    a = 0
    d = -ForwardDiff.gradient(f, start)
    i = 1
    log[:, i] = start

    while i < iterations

        for j = 1:2
            ls(lambda) = f(log[:,i] + lambda * d)

            lambda = newton(ls, newton_start)

            log[:, i+1] = log[:,i] + lambda * d

            a = norm(ForwardDiff.gradient(f, log[:,i+1]))^2 / norm(ForwardDiff.gradient(f, log[:,i]))^2
            d = -ForwardDiff.gradient(f, log[:, i+1]) + a * d
            i = i+1
        end

        d = -ForwardDiff.gradient(f, log[:, i])

        if norm(d) < tolerance
            return log[:,1:i]
        end
    end
    return log
end;

function newton(ls, lambda_0)
    D(f, x) = ForwardDiff.derivative(f, x)
    D2(f, x) = ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), x)

    lambda = lambda_0

    while abs(D(ls, lambda)) > 1e-12
        lambda = lambda - D(ls, lambda) / D2(ls, lambda)
    end

    return lambda
end;

### PROBLEM A ###

f(x) = 0.26 * (x[1]^2 + x[2]^2) - 0.48 * x[1] * x[2]

n = 1000
x = linspace(-10,25,n);
y = linspace(-10,10,n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

contour(x,y,z,
        levels = [0.1, 1, 2, 4 , 7, 10, 15, 20, 30],
        xaxis = (L"$x_1$", (-5,10)),
        yaxis = (L"$x_2$", (-5,10)),
        clims = (0,20),
        clabels = true,
        aspect_ratio = :equal)

gui()

M = 10001 # 10000 iter + start
ϵ = 1e-6 # tolerance
xstart = [7;3] # starting point
n = length(xstart) # dimension of the problem

ts = time();
xj = cjg(f, xstart, M, ϵ, 0)
tend = time() - ts

println("PROBLEM A")
println("Conjugate gradient converged.")
println(" Total time (s): ", tend)
println(" Total steps: ", size(xj, 2)-1)
println(" Sol. found: ", xj[:,end], "/ Opt. value: ", f(xj[:,end]), "\n")

plot!( xj[1,:], xj[2,:], label = "Conjugate", marker=:circle)
savefig("a_path.pdf")

dist_xj = sqrt.(sum(( xj .- [0, 0]).^2,1)');

plot(dist_xj, yscale=:log10, label = "Conjugate",
    xaxis = ("iterations", (1,4)),
    yaxis = (L"$||x_k - \overline{x}||$", ( ϵ, 10)))
savefig("a_convergence.pdf")

###### PROBLEM B #####

f(x) = exp(x[1] + 3*x[2] - 0.1) + exp(x[1] - 3*x[2] - 0.1) + exp(-x[1] - 0.1)

n = 1000
x = linspace(-10,25,n);
y = linspace(-10,10,n);
z = [f([x[i],y[j]]) for j = 1:n, i = 1:n];

contour(x,y,z,
        levels = [2.6, 3, 4, 5 , 6, 7, 8, 10, 12],
        xaxis = (L"$x_1$", (-0.4,-0.3)),
        yaxis = (L"$x_2$", (-0.1,0.1)),
        clims = (0,20),
        clabels = true,
        aspect_ratio = :equal)

gui()


M = 10001 # 10000 iter + start
ϵ = 1e-6 # tolerance
xstart = [1;1.5] # starting point
n = length(xstart) # dimension of the problem

ts = time();
xj = cjg(f, xstart, M, ϵ, 0)
tend = time() - ts;

println("PROBLEM B")
println("Conjugate gradient converged.")
println(" Total time (s): ", tend)
println(" Total steps: ", size(xj, 2)-1)
println(" Sol. found: ", xj[:,end], "/ Opt. value: ", f(xj[:,end]), "\n")

plot!( xj[1,:], xj[2,:], label = "Conjugate", marker=:circle)

savefig("b_path.pdf")

dist_xj = sqrt.(sum(( xj .- [-0.346574, 0]).^2,1)');

plot(dist_xj, yscale=:log10, label = "Conjugate",
    xaxis = ("iterations", (1,10)),
    yaxis = (L"$||x_k - \overline{x}||$", ( ϵ, 10)))

savefig("b_convergence.pdf")
