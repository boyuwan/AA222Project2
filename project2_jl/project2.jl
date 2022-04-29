#=
Working kinda
        AA 222 Project 2, Constrained Optimization
        File Name: project2.jl
        File Description: This file implements two constrained optimization algorithms to optimize the given objective functions. 
        Author: Betty Wan, boyuwan@stanford.edu
        Date: Apr 27th, 2022

        References: TODO
        [1] M. J. Kochenderfer and T. A. Wheeler, Algorithms for Optimization. 
            Textbook contents from Chapter 7.3 Hooke-Jeeves & Algorithm 7.5. The Hooke-Jeeves method were referenced and taken inspiration from.
        [2] AA 222 Project Tutorial, https://www.youtube.com/watch?v=ZPZ9xknXGcQ
=#

#=
    If you want to use packages, please do so up here.
    Note that you may use any packages in the julia standard library
    (i.e. ones that ship with the julia language) as well as Statistics
    (since we use it in the backend already anyway)
=#

# Example:
using LinearAlgebra


#=
    If you're going to include files, please do so up here. Note that they
    must be saved in project1_jl and you must use the relative path
    (not the absolute path) of the file in the include statement.

    [Good]  include("somefile.jl")
    [Bad]   include("/pathto/project1_jl/somefile.jl")
=#

# Example
# include("myfile.jl")


"""
    optimize(f, g, c, x0, n, prob)

Arguments:
    - `f`: Function to be optimized
    - `g`: Gradient function for `f`
    - `c`: Constraint function for 'f'
    - `x0`: (Vector) Initial position to start from
    - `n`: (Int) Number of evaluations allowed. Remember `g` costs twice of `f`
    - `prob`: (String) Name of the problem. So you can use a different strategy for each problem. E.g. "simple1", "secret2", etc.

Returns:
    - The location of the minimum
"""
function optimize(f, g, c, x0, n, prob)
    # if prob == "simple1" || prob == "simple3"
    xhistory, fhistory = penaltyOptimization(f, g, c, x0, n, 0)   # "mixed" = 2, "count" = 1, or "quad" = 0
    x_best = xhistory[argmin(fhistory)]
    return x_best
end


"""
    A helper function that helps construct the i-th basis vector of length n. Referenced from Textbook [1], P99
"""
basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]


function penaltyOptimization(f, g, c, x0, n, type; ρ=30.0, γ=2.5, rho1 = 5.0, rho2 = 10.0, ε = 0.001, ceiling = 10)
    p_quad(x) = sum((max.(0,c(x))).^2)
    p_count(x) = sum(max.(0,c(x)))
    p_mixed(x) = rho1 * p_quad(x) + rho2 * p_count(x)
    f_new(x) = 0

    if type == 0
        f_new(x) = f(x) + ρ*p_quad(x)
    elseif type == 1
        f_new(x) = f(x) + ρ*p_count(x)
    elseif type == 2 
        f_new(x) = f(x) + ρ*p_mixed(x)
    else
        print("Wrong Type! ")
    end

    dim = length(x0)
    xold = zeros(dim,1)
    xnew = zeros(dim,1)
    xold = x0
    xnew = x0 - 0.01.*ones(dim)
    Δ = Inf
    f_newnew = Inf
    gold = g(x0)
    xhistory = [x0]
    fhistory = [f_new(x0)]

    while (count(f, g, c) < (n - 2*dim - 2*dim)) && (sum(abs.(Δ)) > ε) 
        gnew = getGradient(f_new, xnew)
        Δ = ((xnew - xold)/(gnew - gold)) * gnew
        f_newnew = f_new(xnew)
        push!(xhistory, xnew)
        push!(fhistory, f_newnew)
        if f_newnew == f(xnew)
            return xhistory, fhistory 
        end
        if (f_newnew > ceiling)
            xnew = hooke_jeeves(f_new, f, g, c, xnew, n-4*dim, 0.1, 1, 0.5) #(f_new, f, g, c, x, kmax, ε, alpha, γ) 
        end
        xold = xnew
        xnew = xnew -  Δ
        gold = gnew
        ρ *= γ
    end

    return xhistory, fhistory 

end


# Computing gradient with "Central Difference" method
function getGradient(f, x; alpha = 0.0001)
    dim = length(x)
    gradient = zeros(dim)
    for i = 1 : dim
        step = zeros(dim)
        step[i] = alpha
        xpos = x + step
        xneg = x - step
        gradient[i] = (f(xpos) - f(xneg))/(2.0*alpha)
        # print(f(x))
    end
    return gradient
end

function hooke_jeeves(f, f_old, g, c, x, kmax, ε, alpha=1, γ=0.5) 
    y, n = f(x), length(x)
    while (alpha > ε) && (count(f_old, g, c) < kmax-3*n)
        improved = false 
        x_best, y_best = x, y 
        for i in 1 : n
            for sgn in (-1,1)
                x′ = x + sgn*alpha*basis(i, n) 
                y′ = f(x′)
                if y′ < y_best
                    x_best, y_best, improved = x′, y′, true 
                end
            end 
        end
        x, y = x_best, y_best
        if !improved
            alpha *= γ 
        end
    end
    return x 
end

