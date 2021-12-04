#import necessary libraries
include("../algorithms/QR_GS.jl")
using LinearAlgebra
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

### Some matrices to test
#the identity matrix
I_two = Matrix{Float64}(I, 2, 2)
#diagonal matrix, with 2's along diagonal
two_diag = Matrix{Float64}(2I, 2, 2)
#rotation transformation
rot = [cos(90) -sin(90); sin(90) cos(90)]

q1, r1, true_lambdas1, qr_lambdas1, iters1 = qr_eigenvals(I_two)
my_display(q1, r1, true_lambdas1, qr_lambdas1, iters1, "I_two")

q2, r2, true_lambdas2, qr_lambdas2, iters2 = qr_eigenvals(two_diag)
my_display(q2, r2, true_lambdas2, qr_lambdas2, iters2, "two_diag")

q3, r3, true_lambdas3, qr_lambdas3, iters3 = qr_eigenvals(rot)
my_display(q3, r3, true_lambdas3, qr_lambdas3, iters3, "rot")