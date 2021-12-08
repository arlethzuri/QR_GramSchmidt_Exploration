"""
Identity Matrices:

Here, I wanted to confirm for myself that the QR algorithm takes the same number of iterations
to determine the eigenvalues for any size identity matrix.

A plot of iterations QR runs as identity matrix size changes is generated.
"""
#import necessary libraries
include("../../algorithms/QR_GS.jl")
using LinearAlgebra
using Plots
using .QR_GS

#display items for each matrix
function my_display(q, r, true_lambdas, qr_lambdas, iters, display_name)
    println("------------------------")
    println("For the matrix ", display_name)
    println("Q is: ")
    display(q)
    println("")
    
    println("R is: ")
    display(r)
    println("")
    
    println("true eigenvalues are: ")
    display(true_lambdas)
    println("")
    
    println("Eigenvalues determined by QR algorithm are: ")
    display(qr_lambdas)
    println("")

    println("Iterations: ", iters)
    println("------------------------")
end
#max identity matrix size to test to
max_size = 100

#iteration count holder for each matrix size, I go up to 100
iterations = zeros(max_size)

#create identity matrices of sizes 1 to 100 and run through QR algorithm
for i = 1:max_size
    #create identity matrix
    curr_matrix = Matrix(I, i, i)

    #run through QR algorithm, get Q, R, "true" eigvals, 
    # QR estimated eigvals, and iter count
    q, r, true_lambdas, qr_lambdas, iters = qr_eigenvals(curr_matrix)

    #store iteration count in iterations vector
    iterations[i] = iters
end

plot(iterations,
    xlabel = "Size N of N×N Identity Matrix",
    ylabel = "Iteration Count",
    label = "Iterations",
    title = "QR Iterations for Growing Identity Matrix Size",
    fontfamily = "Times",
    dpi=1000)
savefig("plots/iterations")