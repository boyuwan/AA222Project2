#=
        AA 222 Project 1, Unconstrained Optimization
        File Name: plotting.jl
        File Description: This file contains the implementation of plots required for Assignment 1 submission.
        Author: Betty Wan, boyuwan@stanford.edu
        Date: Apr 11st, 2022

        References: 
        [1] M. J. Kochenderfer and T. A. Wheeler, Algorithms for Optimization. 
            Textbook contents from Chapter 7.3 Hooke-Jeeves & Algorithm 7.5. The Hooke-Jeeves method were referenced and taken inspiration from.
        [2] AA 222 Project Tutorial, https://www.youtube.com/watch?v=ZPZ9xknXGcQ
=#

# import Pkg
# Pkg.add("Plots")
using Plots

include("project2_jl/helpers.jl")
include("project2_jl/simple.jl")
include("project2_jl/project2.jl")

x_rose_1, f_rose_1 = penalty_method(simple1, simple1_gradient, simple1_constraints, [1.4,0.4], 2000, "count") 
plot(collect(1:length(f_rose_1)), f_rose_1, xlabel = "Iteration", ylabel = "f(x)",  label="simple1, x0 = [1.4, 0.4]")
savefig("Simple 1 Convergence.png")
print(f_rose_1)

# basis(i, n) = [k == i ? 1.0 : 0.0 for k in 1 : n]


# function rosenbrock(x,y)
#     return (1.0 - x)^2 + 100.0 * (y - x^2)^2
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
#                 # if count(f, g) >= n
#                 #     return xhistory, fhistory
#                 # end
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

# # Plotting optimizer function
# x_rose_1, f_rose_1 = hookeJeeves(rosenbrock, rosenbrock_gradient, [1.4,0.4], 100, 0.10, 0.001, 0.5)  
# x_rose_2, f_rose_2 = hookeJeeves(rosenbrock, rosenbrock_gradient, [1.8,1.2], 100, 0.10, 0.001, 0.5)  
# x_rose_3, f_rose_3 = hookeJeeves(rosenbrock, rosenbrock_gradient, [0.4,1.2], 100, 0.10, 0.001, 0.5) 

# x_himm_1, f_himm_1 = hookeJeeves(himmelblau, himmelblau_gradient, [1.4,0.4], 100, 0.10, 0.001, 0.5)  
# x_himm_2, f_himm_2 = hookeJeeves(himmelblau, himmelblau_gradient, [1.8,1.2], 100, 0.10, 0.001, 0.5)  
# x_himm_3, f_himm_3 = hookeJeeves(himmelblau, himmelblau_gradient, [2.4,1.2], 100, 0.10, 0.001, 0.5) 

# x_pow_1, f_pow_1 = hookeJeeves(powell, powell_gradient, [1.4, 0.4, 1.5, 1.5], 100, 0.10, 0.001, 0.5)  
# x_pow_2, f_pow_2 = hookeJeeves(powell, powell_gradient, [1.8, 1.2, 2, 2], 100, 0.10, 0.001, 0.5)  
# x_pow_3, f_pow_3 = hookeJeeves(powell, powell_gradient, [2.4, 1.2, 1, 1], 100, 0.10, 0.001, 0.5) 

# # Contour Plot of Rosenbrock
# xr = -2:0.1:2
# yr = -2:0.1:2
# # println(x_rose_3)

# contour(xr, yr, rosenbrock, levels = [10,25,50,100,200,250,300], colorbar = false, c = cgrad(:viridis, rev = true), legend = false, xlims = (-2, 2), ylims = (-2, 2), xlabel = "x1", ylabel = "x2", aspectratio = :equal, clim = (2, 500))
# plot!([x_rose_1[i][1] for i = 1:length(x_rose_1)], [x_rose_1[i][2] for i = 1:length(x_rose_1)], color = :black, label="Hooke Jeeves Method, x0 = [1.4, 0.4]")
# plot!([x_rose_2[i][1] for i = 1:length(x_rose_2)], [x_rose_2[i][2] for i = 1:length(x_rose_2)], color = :blue, label="Hooke Jeeves Method, x0 = [1.8, 1.2]")
# plot!([x_rose_3[i][1] for i = 1:length(x_rose_3)], [x_rose_3[i][2] for i = 1:length(x_rose_3)], color = :red, label="Hooke Jeeves Method, x0 = [0.4, 1.2]")
# savefig("Rosenbrock Contour Plot and Optimization Paths.png")

# # Convergence Plot of Rosenbrock
# plot(collect(1:length(f_rose_1)), f_rose_1, xlabel = "Iteration", ylabel = "f(x)",  label="Hooke Jeeves Method, x0 = [1.4, 0.4]")
# plot!(collect(1:length(f_rose_2)), f_rose_2, label="Hooke Jeeves Method, x0 = [1.8, 1.2]")
# plot!(collect(1:length(f_rose_3)), f_rose_3, label="Hooke Jeeves Method, x0 = [0.4, 1.2]")
# savefig("Rosenbrock Optimization Convergence.png")

# # Convergence Plot of Himmelblau
# plot(collect(1:length(f_himm_1)), f_himm_1, xlabel = "Iteration", ylabel = "f(x)",  label="Hooke Jeeves Method, x0 = [1.4, 0.4]")
# plot!(collect(1:length(f_himm_2)), f_himm_2, label="Hooke Jeeves Method, x0 = [1.8, 1.2]")
# plot!(collect(1:length(f_himm_3)), f_himm_3, label="Hooke Jeeves Method, x0 = [2.4, 1.2]")
# savefig("Himmelblau Optimization Convergence.png")

# # Convergence Plot of Powell
# plot(collect(1:length(f_pow_1)), f_pow_1, xlabel = "Iteration", ylabel = "f(x)",  label="Hooke Jeeves Method, x0 = [1.4, 0.4, 1.5, 1.5]")
# plot!(collect(1:length(f_pow_2)), f_pow_2, label="Hooke Jeeves Method, x0 = [1.8, 1.2, 2, 2]")
# plot!(collect(1:length(f_pow_3)), f_pow_3, label="Hooke Jeeves Method, x0 = [2.4, 1.2, 1, 1]")
# savefig("Powell Optimization Convergence.png")