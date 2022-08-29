
from sys import exit
import sys
import math
import spkmeans_module as spkm
import pandas as pd
import numpy as np


def invalid_input():
    print("Invalid Input!")
    exit()


def other_error():
    print("An Error Has Occurred")
    exit()


def is_valid_path(path):
    """Function to check if the given path is valid - ending with .txt or .csv"""
    if len(path.split('.')) == 2:
        if path.split('.')[1] == 'txt' or path.split('.')[1] == 'csv':
            return True
    return False


def print_2d_list(arr):
    """Simple function to print 2d array"""
    for row in arr:
        for i in range(len(row)):
            row[i] = f'{row[i]:.4f}'
        print(','.join(row))


def k_means_pp(points, K):
    """Function whiche calculate and return the list of inital centroids for Kmeans algorithm"""
    centroids = np.zeros((K, points.shape[1]))
    centroids_indexes = [0] * K
    points_number = points.shape[0]
    D_values = np.zeros(points_number)
    P_values = np.zeros(points_number)
    centroids_indexes[0] = np.random.choice(len(points))
    centroids[0] = points[centroids_indexes[0]]
    for i in range(1, K):
        for l in range(points_number):

            distance_list = [
                (points[l] - centroids[j]).dot((points[l] - centroids[j])) for j in range(i)]
            D_values[l] = min(distance_list)
        for l in range(points_number):
            P_values[l] = D_values[l] / sum(D_values)
        centroids_indexes[i] = np.random.choice(len(points), p=P_values)
        centroids[i] = points[centroids_indexes[i]]
    print(','.join([f'{c}' for c in centroids_indexes]))
    return centroids


def spk(points, K):
    """Spectral clustering algorithm. it prints the centroids"""
    # step 1: calculate wam(points)
    diagonal = spkm.jacobi_lnorm(points)
    narr = np.array(diagonal)
    df = pd.DataFrame(narr.T)
    df = df.sort_values(by=[0], ascending=False).reset_index(drop=True)

    if K == 0:
        # determine K
        delta_max = 0
        for i in range(0, math.floor(len(diagonal[0]) / 2)):
            delta = abs(df[0][i] - df[0][i + 1])
            if delta > delta_max:
                K = i + 1
                delta_max = delta

    # construct U - k largest eiganvectors of Lnorm in ascedning order
    U = df.iloc[:K, :].iloc[:, 1:].to_numpy().T
    # T is normalized U
    for i in range(U.shape[0]):
        norm = np.linalg.norm(U[i])
        if norm:  # checking if the norm isn't 0
            U[i] /= norm
        T = U
    # calling kmeans algorithm with T
    centroids = spkm.kmeans_fit(K, 300, T.tolist(), k_means_pp(T, K).tolist(), T.shape[0], T.shape[1], 0)
    print_2d_list(centroids)
    return(centroids)


def main():

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

    if(input_goal == "spk"):
        spk(points.tolist(), K)
    elif(input_goal == "wam"):
        print_2d_list(spkm.wam_fit(points.tolist()))
    elif(input_goal == "ddg"):
        print_2d_list(spkm.ddg_fit(points.tolist()))
    elif(input_goal == "lnorm"):
        print_2d_list(spkm.lnorm_fit(points.tolist()))
    elif(input_goal == "jacobi"):
        print_2d_list(spkm.jacobi_fit(points.tolist()))
    else:
        invalid_input()


if __name__ == '__main__':
    main()
