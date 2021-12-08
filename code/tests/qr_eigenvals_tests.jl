#import necessary libraries
include("../algorithms/QR_GS.jl")
using LinearAlgebra
using .QR_GS

### Some matrices to test
#the identity matrix
I_two = Matrix{Float64}(I, 2, 2)
#diagonal matrix, with 2's along diagonal
two_diag = Matrix{Float64}(2I, 2, 2)
#rotation transformation
rot = [cos(90) -sin(90); sin(90) cos(90)]

q1, r1, true_lambdas1, qr_lambdas1, iters1, errors1 = qr_eigenvals(I_two)
my_display(q1, r1, true_lambdas1, qr_lambdas1, iters1, errors1, "I_two")

q2, r2, true_lambdas2, qr_lambdas2, iters2, errors2 = qr_eigenvals(two_diag)
my_display(q2, r2, true_lambdas2, qr_lambdas2, iters2, errors2, "two_diag")

q3, r3, true_lambdas3, qr_lambdas3, iters3, errors3 = qr_eigenvals(rot, 1e-10, 10)
my_display(q3, r3, true_lambdas3, qr_lambdas3, iters3, errors3, "rot")