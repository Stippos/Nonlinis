@enum LS NEWTON ARMIJO

using ForwardDiff



function Newton(g, x0, N, LS, ϵ)
    counter = 0;
    p = x0;
    for k = 1:N
        gd = -ForwardDiff.hessian(g, p) \ ForwardDiff.gradient(g, p);
        if norm(gd) < ϵ
            break
        else
            ls_function(lambda) = g(p + lambda * gd);
            if LS == NEWTON
                lambda = newton_ls(ls_function, ϵ);
            elseif LS == ARMIJO
                lambda = armijo(ls_function, p, gd);
            end;
        end;
        p = p + lambda * gd;
        counter += 1;
    end;
    return g(p), counter;
end;

function BFGS(g, x0, N, LS, ϵ)
    p = x0;
    B = ForwardDiff.hessian(g, p);
    counter = 0;

    for k = 1:N
        #println(B)
        gd = -inv(B) * ForwardDiff.gradient(g, p)
        #println(gd)
        ls_function(lambda) = g(p + lambda * gd);

        if LS == NEWTON
            lambda = newton_ls(ls_function, ϵ);
        elseif LS == ARMIJO
            lambda = armijo(ls_function, p, gd);
        end;

        step = lambda * gd
        p = p + step;
        counter += 1;

        if norm(ForwardDiff.gradient(g, p)) < ϵ
            break
        end

        q = ForwardDiff.gradient(g, p) - gd;
        #println(lambda);
        ##println(gd);
        #println(step);
        #println(q);
        B = B + (q*q') / (q'*step) - (B*step*step'*B) / (step'*B*step);

    end;
    return g(p), counter ;
end

function Gradient(g, x0, N, LS, ϵ)
    counter = 0;
    p = x0;

    for k = 1:N
        gd = -ForwardDiff.gradient(g, p);
        if norm(gd) < ϵ
            break
        else
        ls_function(lambda) = g(p + lambda * gd);

            if LS == NEWTON
                lambda = newton_ls(ls_function, ϵ);
            elseif LS == ARMIJO
                lambda = armijo(ls_function, p, gd);
            end;

            p = p + lambda * gd;
            counter += 1;
        end;
    end;

    return g(p), counter;

end;

function newton_ls(h, ϵ)

    derivative(h, lambda) = ForwardDiff.derivative(h, lambda);
    second_derivative(h, lambda) = ForwardDiff.derivative(y -> ForwardDiff.derivative(h, y), lambda);

    lambda = 1;

    while abs(derivative(h, lambda)) > 1e-10
        next_lambda = lambda - derivative(h, lambda)/ second_derivative(h, lambda);
        lambda = next_lambda;
    end;
    return lambda;
end;

function armijo(fun, x, d)
    alpha = 0.1;
    beta = 0.7;
    lambda = 10;
    slope = ForwardDiff.derivative(fun, 0);
    x0 = fun(0);

    while fun(lambda) > x0 + alpha * lambda * slope
        lambda = beta * lambda;
    end;

    return lambda;

end;
#Pkg.add("ForwardDiff")
