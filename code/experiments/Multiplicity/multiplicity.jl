"""
Noise on a Nice Diagonal Matrix

In this experiment, I want to see how the QR algorithm performs on a diagonal matrix that progressively has
more noise added to it. By noise, I mean going from a nice diagonal matrix where all entries are 0 except the diagonal
to a diagonal matrix that has noise added to the entire matrix as demonstrated below.

A plot showing iteration count by noise addition to matrix is generated.

"""
#import necessary libraries
include("../../algorithms/QR_GS.jl")
using LinearAlgebra
using StableRNGs
using Plots
using .QR_GS

#seed for random numbers
rng = StableRNGs.StableRNG(1)

#iterations to perform
iters_todo = 2
nice_size = 15

#to hold iteration count
iterations = zeros(iters_todo)

#create the nice diagonal matrix of size 100,
# simply k = 1, 2, ... 100 along diagonal of 0 matrix
nice = diagm(0 => 1:1:nice_size)

#run this nice matrix through QR algorithm
q1, r1, true_lambdas1, qr_lambdas1, iters1, _ = qr_eigenvals(nice)
iterations[1] = iters1

#progressively add more noise, 10 times
for i = 2:iters_todo
    println("On iteration ", i)
    #generate noise
    noise = i .* rand(rng, nice_size, nice_size)
    #add noise to nice
    noisy = nice + noise
    #run matrix through QR algorithm
    q, r, true_lambdas, qr_lambdas, iters, errors = qr_eigenvals(noisy, 1e-5, 500)
    display(log.(errors))
    #plot for noise level i
    plot(log.(errors),
        xlabel = "Iterations",
        ylabel = "Error between QR and Julia Eigenvalues",
        label = "L2 Norm Difference",
        title = string("QR Error for Noise Level ", i, "× 1e-4"),
        fontfamily = "Times",
        dpi=1000)
    savefig(string("plots/noise_lvl_", i))

    iterations[i] = iters
end

#iteration plot
plot(iterations,
    xlabel = "Noise of magnitude (1e-5 x N) added",
    ylabel = "Iteration Count",
    label = "Iterations",
    title = "QR Iterations for Increasing Noise",
    fontfamily = "Times",
    dpi=1000)
savefig("plots/noise")