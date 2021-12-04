module QR_GS

export qr_gs, qr_eigenvals

include("QR_GSsource.jl")

"""
Compute an orthonormal basis for a matrix of any size,
Inputs:
    x - matrix to find orthonormal basis of
Outputs:
    u - orthonormal basis of x
"""
function qr_gs(x)
    qr_gs_(x)
end

"""
Run the QR algorithm on a full-rank, real, square matrix x and iterate until
 eigenvalues within tolerance level are found or iteration limit has been reached
Inputs:
    x - matrix to find eigenvalues of using QR
    tol - tolerance level
    max_iters - iteration limit
Outputs:
    q - the Q matrix of x's QR decomposition
    r - the R matrix of x's QR decomposition
    true_lambdas - the true eigenvalues of x according to Julia
    qr_lambdas - the eigenvalues determined by QR algorithm
    iters - the number of iterations x ran through the QR algorithm
"""
function qr_eigenvals(x, tol::Float64=default_tol, max_iters::Int64=default_maxit)
    qr_eigenvals_(x, tol, max_iters)
end

end #module