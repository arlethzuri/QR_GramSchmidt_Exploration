#import necessary libraries
using LinearAlgebra

default_maxit = 1000000
default_tol = 1e-10

"""
Determine if a given matrix is full-rank and square
Inputs:
    x - matrix to check
Outputs:
    frrs - boolean saying whether x is full-rank, real, and square
"""
function fullrank_square(x)
    #get size of matrix
    rows, cols = size(x)

    #if rank of given matrix does not equal dim(Col(X)), throw error
    if rank(x) != cols
        return false
    end

    #if matrix given is not square throw error
    if rows != cols
        return false
    end

    return true
end

"""
Compute an orthonormal basis for a full-rank, real, square matrix of any size
Inputs:
    x - matrix to find orthonormal basis and QR decomposition of
Outputs:
    q - Q matrix in QR decomposition
    r - R matrix in QR decomposition
"""
function qr_gs_(x)
    #check if matrix meets pre-requisites
    if fullrank_square(x) == false
        throw("Matrix is not square or it is not full-rank!")
        return
    end

    #get size of matrix
    rows, cols = size(x)

    #initialize orthogonal matrix and upper triangular matrix
    q = zeros(rows, cols)
    r = zeros(rows, cols)

    #define v1 & e1, this is the first step before we can begin iterating
    v1 = x[:, 1]
    e1 = v1 ./ norm(v1)

    #set the first columns of q and r
    q[:, 1] = e1
    a_i = dot(e1, x[:, 1]) #variable to hold last sum a_i, we begin with a_1
    r[1, 1] = a_i

    #continue to iterate through gram-schmidt through all columns, fill q & r
    for i = 2:cols
        #get the ith column from the x matrix and e_i
        x_i = x[:, i]
        #initialize weighted sum of projections
        projection_sum = zeros(rows)

        #calculate the v_i-th column of the orthonormal basis
        for j = 1:i-1
            curr_q = q[:, j]
            projection_sum += (dot(curr_q, x_i)/norm(curr_q)) .* curr_q
        end
        #subtract projections and normalize
        e_i = (x_i - projection_sum)
        e_i ./= norm(e_i)

        q[:, i] = e_i

        #fill in the r matrix with r_entry
        for k = 1:rows
            r[k, i] = dot(x_i, q[:, k])
        end
    end

    #return the orthonormal basis, q, and r
    return q, r
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
function qr_eigenvals_(x, tol::Float64=default_tol, max_iters::Int64=default_maxit)
    #check if matrix meets pre-requisites
    if fullrank_square(x) == false
        throw("Matrix is not square or it is not full-rank!")
        return
    end

    #get size of matrix
    rows, cols = size(x)
    
    #use Julia's eigenvalue solver to get what we will consider the "true"
    # eigenvalues of x
    true_lambdas = eigvals(x)

    #initialize vector to hold lambda estimates
    qr_lambdas = zeros(cols)

    #initialize matrices q and r
    q = zeros(rows, cols)
    r = zeros(rows, cols)

    #current error, initialize iteration counter
    err = 1
    iters = 0

    #initialize matrix to find q r decomposition of to be original matrix x
    curr_A = x

    #find QR decomposition of x and iterate until max_iters have occured
    # or tolerance level has been satisfied
    while (err > tol) && (iters < max_iters)
        q, r = qr_gs_(curr_A)
        curr_a = r * q

        #get values along diagonal of curr_A
        qr_lambdas = diag(curr_A)
        #calculate error, using L2 norm
        err = norm(true_lambdas) - norm(qr_lambdas)
        #increment iteration counter
        iters += 1
    end

    #return the final qr decomposition, "true" eigenvalues, estimated eigenvalues and iteration count
    return q, r, true_lambdas, qr_lambdas, iters
end