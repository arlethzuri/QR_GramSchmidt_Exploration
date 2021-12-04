#import necessary libraries
include("../algorithms/QR_GS.jl")
using LinearAlgebra
using .QR_GS

#using example from https://rpubs.com/aaronsc32/qr-decomposition-gram-schmidt
A = [2 -2 18; 2 1 0; 1 2 0]
B = [0 3; 4 1; 0 4]

Q, R = qr_gs(A)

println("Q is: ")
display(Q)
println("")
println("R is: ")
display(R)
println("")