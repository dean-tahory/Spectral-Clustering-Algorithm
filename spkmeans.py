
import sys
import math
import spkmeans as spkmeans
import pandas as pd
import numpy as np


def invalid_input():
    print("Invalid Input!")
    exit()


def other_error():
    print("An Error Has Occurred")
    exit()


def is_valid_path(path):
    if len(path.split('.')) == 2:
        if path.split('.')[1] == 'txt':
            return True
    return False


def print_2d_list(arr):
    for row in arr:
        for i in range(len(row)):
            row[i] = f'{row[i]:.4f}'
        print(','.join(row))


# first step: validate command line foramt :
if len(sys.argv) == 4:
    if not sys.argv[1].isnumeric():
        invalid_input()
    K = int(sys.argv[1])
    input_goal = sys.argv[2]
    input_file = sys.argv[3]

    if not is_valid_path(input_file):
        invalid_input()
else:
    invalid_input()

# second step: reading the input, checking if K < N
try:
    points_table = pd.read_csv(input_file, header=None, dtype=float)

except FileNotFoundError:
    other_error()
if(K >= len(points_table.index)):
    invalid_input()
points = points_table.to_numpy()
np.random.seed(0)


# third step: k-means++ algorithm
max_iter = 300
eps = 0.001


def k_means_pp(points, K):
    centroids = np.zeros((K, points.shape[1]))
    centroids_indexes = [0]*K
    points_number = points.shape[0]
    D_values = np.zeros(points_number)
    P_values = np.zeros(points_number)
    centroids_indexes[0] = np.random.choice(len(points))
    centroids[0] = points[centroids_indexes[0]]
    for i in range(1, K):
        for l in range(points_number):

            distance_list = [
                (points[l]-centroids[j]).dot((points[l]-centroids[j])) for j in range(i)]
            D_values[l] = min(distance_list)
        for l in range(points_number):
            P_values[l] = D_values[l]/sum(D_values)
        centroids_indexes[i] = np.random.choice(len(points), p=P_values)
        centroids[i] = points[centroids_indexes[i]]
    print(','.join([f'{c}' for c in centroids_indexes]))
    return centroids


def determine_k(matrix, dim):
    # finding the eiganvalues of the matrix with jacobi method
    eiganvalues = spkmeans.jacobi_fit(matrix, dim)[0]
    eiganvalues.sort(reverse=True)
    k = 0
    delta_max = 0
    for i in range(0, math.floor(dim/2)):
        delta = abs(eiganvalues[i]-eiganvalues[i+1])
        if delta > delta_max:
            k = i
            delta_max = delta
    return k+1

# TEST


def spk(points: list[float], dim: float, K: int):
    # step 1: calculate wam(points)
    lnorm_matrix = spkmeans.lnorm_fit(points, dim)
    if K == 0:
        K = determine_k(lnorm_matrix, dim)
    eigvalue_to_eigvector = {}

    # TODO change points to lnorm_matrix
    jacobi_mat = spkmeans.jacobi_fit(lnorm_matrix, dim)
    for i in range(dim):
        eigvalue_to_eigvector[jacobi_mat[0][i]] = np.array(
            [jacobi_mat[j][i] for j in range(1, dim+1)])
    eiganvalues = jacobi_mat[0]
    eiganvalues.sort()
    U_matrix = eigvalue_to_eigvector[eiganvalues[0]]
    for i in range(1, K):
        U_matrix = np.vstack([U_matrix, eigvalue_to_eigvector[eiganvalues[i]]])
    U_matrix = np.transpose(U_matrix)
    # print_2d_list(U_matrix.tolist())
    V_matrix = U_matrix
    # TEST if V is normalized by rows
    norm_vector = np.array([np.linalg.norm(U_matrix[i])
                           for i in range(U_matrix.shape[0])])
    for i in range(U_matrix.shape[0]):
        V_matrix[i] = U_matrix[i]/norm_vector[i]
    # print('\n')
    # print_2d_list(V_matrix.tolist())
    centroids = spkmeans.kmeans_fit(K, max_iter, V_matrix.tolist(), k_means_pp(
        V_matrix, K).tolist(), V_matrix.shape[0], V_matrix.shape[1], eps)
    for c in centroids:
        for i in range(len(c)):
            c[i] = f'{c[i]:.4f}'
        print(','.join(c))


    # main methods
if(input_goal == "spk"):
    spk(points.tolist(), points.shape[0], K)
elif(input_goal == "wam"):
    print_2d_list(spkmeans.wam_fit(points.tolist(),
                                   points.shape[0]))
elif(input_goal == "ddg"):
    print_2d_list(spkmeans.ddg_fit(points.tolist(),
                  points.shape[0]))
elif(input_goal == "lnorm"):
    print_2d_list(spkmeans.lnorm_fit(points.tolist(),
                  points.shape[0]))
elif(input_goal == "jacobi"):
    print_2d_list(spkmeans.jacobi_fit(points.tolist(), points.shape[0]))
else:
    invalid_input()
