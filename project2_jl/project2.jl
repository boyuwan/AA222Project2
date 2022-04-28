#=
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
    # if prob == "simple1" 
        xhistory, fhistory = penalty_method(f, g, c, x0, n,"quadratic") 
    # elseif prob == "simple2"
    #     xhistory, fhistory = hookeJeeves(f, g, x0, n, 1.5, 0.001, 0.5)
    # elseif prob == "simple3"
    #     xhistory, fhistory = hookeJeeves(f, g, x0, n, 0.5, 0.001, 0.5)
    # else 
    #     xhistory, fhistory = hookeJeeves(f, g, x0, n, 1, 0.001, 0.5)
    # end
    # print("hi")
    x_best = xhistory[argmin(fhistory)]
    return x_best
end


"""
    A helper function that helps construct the i-th basis vector of length n. Referenced from Textbook [1], P99
"""
basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]

function pcount(c,x,p)
    constraints = c(x)
    dim = length(constraints)
    for i in 1 : dim
        violated = constraints[i] > 0 ? 1 : 0
        p += violated 
    end
    return p
end

function pquad(c,x,p)
    constraints = c(x)
    dim = length(constraints)
    for i in 1 : dim
        violated = constraints[i] > 0 ? constraints[i] : 0
        p += violated * violated
    end
    return p
end


function minimize(f,x, state)
    if state == "hj"
        return hooke_jeeves(f, x) 
    end
end

function penalty_method(f, g, c, x, k_max,type; ρ=1, γ=2) 
    
    # Initialization
    p = 0
    xhistory = [x]
    fhistory = [f(x)]

    if type == "count"
        while (count(f, g, c) < k_max-1)
            fnew = f(x)
            p = pcount(c,x,p)
            x = minimize(x -> (fnew + ρ*p), x, "hj") 
            ρ *= γ
            push!(xhistory, x)
            push!(fhistory, fnew)
            if p == 0
                return xhistory, fhistory 
            end
        end
    elseif type == "quadratic"
        while (count(f, g, c) < k_max-1)
            fnew = f(x)
            p = pquad(c,x,p)
            x = minimize(x -> fnew + ρ*p, x, "hj") 
            ρ *= γ
            push!(xhistory, x)
            push!(fhistory, fnew)
            if p == 0
                return xhistory, fhistory 
            end
        end
    end

    # Reset values
    p = 0
    return xhistory, fhistory 

end

function hooke_jeeves(f, x; alpha=1, ε=0.0001, γ=0.5) 
    y, n = f(x), length(x)
    while alpha > ε
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


# function getP(g, x)
#     dim = length(x)
#     p = 0
#     gfunc = g(x) 
#     for i in 1 : dim
#         p -= (1/gfunc[i])
#     end
#     return p
# end

"""
    interior_point_method(f, g, x; ρ=1, γ=2, ε=0.001)    

    An optimization function that performs direct method optimization on a given objective function. Referenced from Textbook [1], P104
Arguments:
    - `f`: Function to be optimized
    - `g`: Gradient function for `f`
    - `x`: (Vector) Initial position to start from
    - `n`: (Int) Number of evaluations allowed.
    - `ρ`: The initial penalty, ρ > 0
    - `ε`: The stopping tolerance, ε > 0
    - `γ`: Penalty Multiplier, γ > 1

Returns:
    - The history of x and f explored
"""
# function interior_point_method(f, p, x, n; ρ=1, γ=2, ε=0.001) 

#     minimizeQuadratic
#     delta = Inf
#     xhistory = [x]
#     fhistory = [f(x)]

#     while (delta > ε) && (count(f, g, c) < n)
#         pnew = getP(g,x)
#         xnew = minimize(x -> f(x) + pnew/ρ, x) 
#         delta = norm(xnew - x)
#         x = xnew
#         push!(xhistory, xnew)
#         push!(fhistory, ynew)
#         ρ *= γ
#     end
#     return xhistory, fhistory  
# end

# function hookeJeeves(f, g, x, n, a, ε, γ=0.5) 

#     xhistory = [x]
#     fhistory = [f(x)]
#     y, dim = fhistory[1], length(x)

#     while a > ε

#         improved = false 
#         x_best, y_best = x, y

#         for i in 1 : dim
#             for sgn in (-1,1)
#                 if count(f, g) >= n
#                     return xhistory, fhistory
#                 end
#                 xnew = x + sgn * a * basis(i, dim) 
#                 ynew = f(xnew)
#                 if ynew < y_best
#                     push!(xhistory, xnew)
#                     push!(fhistory, ynew)
#                     x_best, y_best, improved = xnew, ynew, true
#                 end
#             end 
#         end
#         x, y = x_best, y_best

#         if !improved
#             a *= γ 
#         end
    
#     end
#     return xhistory, fhistory 

# end