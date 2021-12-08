"""
Imaginary Matrices:

Here, I wanted to see how the QR algorithm would react to working on matrices with strictly imaginary
eigenvalues.

Plots de
"""
#import necessary libraries
include("../../algorithms/QR_GS.jl")
using LinearAlgebra
using Plots
using .QR_GS

##Some matrices to try
#2x2 strictly imaginary eigenvalued matrix, the rotation matrix
rot = [cos(90) -sin(90); sin(90) cos(90)]
m1 = [1 -1; 1 1]
m2 = [4/5 -3/5 0; 3/5 4/5 0; 1 2 2]
m3 = [-1 -4; 1 -1]
m4 = [-3 5; -2 3]

q1, r1, true_lambdas1, qr_lambdas1, iters1, errors1 = qr_eigenvals(rot, 1e-10, 10)
my_display(q1, r1, true_lambdas1, qr_lambdas1, iters1, errors1, "rot")

q2, r2, true_lambdas2, qr_lambdas2, iters2, errors2 = qr_eigenvals(m1, 1e-10, 10)
my_display(q2, r2, true_lambdas2, qr_lambdas2, iters2, errors2, "m1")

q3, r3, true_lambdas3, qr_lambdas3, iters3, errors3 = qr_eigenvals(m2, 1e-10, 10)
my_display(q3, r3, true_lambdas3, qr_lambdas3, iters3, errors3, "m2")

q4, r4, true_lambdas4, qr_lambdas4, iters4, errors4 = qr_eigenvals(m3, 1e-10, 10)
my_display(q4, r4, true_lambdas4, qr_lambdas4, iters4, errors4, "m3")

q5, r5, true_lambdas5, qr_lambdas5, iters5, errors5 = qr_eigenvals(m4, 1e-10, 10)
my_display(q5, r5, true_lambdas5, qr_lambdas5, iters5, errors5, "m4")