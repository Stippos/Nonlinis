@enum LS NEWTON ARMIJO

using ForwardDiff

function Newton(g, x0, N, LS, ϵ)
    return -1, -1
end

function BFGS(g, x0, N, LS, ϵ)
    return -1, -1
end

function Gradient(g, x0, N, LS, ϵ)
    iters = 1000
    counter = 0
    point = x0
    gd = -ForwardDiff.gradient(g, point)

    for k = 1:iters
        if norm(gd) < ϵ
            break
        else



    return g(point), counter

    return -1, -1
end

Pkg.add("ForwardDiff")
