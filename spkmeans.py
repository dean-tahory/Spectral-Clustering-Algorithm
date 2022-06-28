
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

# correct input format: K max_iter(optional) eps input_file1.txt input_file2.txt


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
    points_table = pd.read_csv(input_file, header=None)
except FileNotFoundError:
    other_error()

if(K >= len(points_table.index)):
    invalid_input()
points = points_table.to_numpy()
np.random.seed(0)


# main methods
def wam():
    wam_fit()


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
    print(' '.join([f'{c}' for c in centroids_indexes]))
    return centroids


# fourth step: writing output
centroids = spkmeans.fit(K, max_iter, points.tolist(), k_means_pp(
    points, K).tolist(), points.shape[0], points.shape[1], eps)
for c in centroids:
    for i in range(len(c)):
        c[i] = f'{c[i]:.4f}'
    print(','.join(c))
