'''
Gram-Schmidt Orthogalization
'''
import numpy as np

def gram_schmidt(x):
    num_rows = x.shape[0]
    num_cols = x.shape[1]
    cols_to_delete = [] #hold which columns are 0s
    zeros = np.zeros(num_rows, dtype=int) #define a column that has 0s

    #check each column of x
    for i in range(num_cols):
        col = x[:, i] #get ith column
        #if the ith column is all zeros, add i to list
        if np.array_equal(col, zeros):
            cols_to_delete.append(i)

    #delete columns that're all 0s
    x_no_zeros = np.delete(x, cols_to_delete, axis=1)
    #get u1 using the first column of x_no_zeros
    x1 = x_no_zeros[:, 0]
    u = x1/np.linalg.norm(x1, 2)
    #promote u to have another axis
    u = u[:, np.newaxis]

    #if there is more than one column in x
    if x_no_zeros.shape[1] > 1:
        #find the rest of the columns of u
        for j in range(1, x_no_zeros.shape[1]):
            xj_weighted_sum = 0
            print("lmao: ", xj_weighted_sum)
            #find the best rep. of xj as a weighted sum of u_1...u_j-1
            for i in range(j):
                xj_weighted_sum += np.dot(np.dot(u[:, i].T, x_no_zeros[:, j]), u[:, i])
                print("------")
                print(xj_weighted_sum)
                print("------")

            #calculate xj' and normalize it to get jth column of u and append to u
            xj_prime = x_no_zeros[:, j] - xj_weighted_sum
            uj = xj_prime/np.linalg.norm(xj_prime, 2)

            #promote uj to have another axis
            uj = uj[:, np.newaxis]
            #add uj as the jth column of u
            u = np.append(u, uj, axis=1)

    #delete zero columns from u
    rows_to_delete = []
    num_rows = u.shape[0]
    for i in range(num_rows):
        row = u[i] #get ith column
        #if the ith column is all zeros, add i to list
        if np.array_equal(row, zeros):
            rows_to_delete.append(i)

    #delete columns that're all 0s
    u = np.delete(u, rows_to_delete, axis=0)
    
    #return u
    return u

def main():
    #an example from lecture 7
    x = np.array([[1, 1], [1, 3], [0, 0]])
    u = gram_schmidt(x)
    print(u)

main()