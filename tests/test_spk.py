

from sklearn.cluster import SpectralClustering
import matplotlib.pyplot as plt
import pandas as pd


import sys
sys.path.append("../Spectral-Clustering-Algorithm/src")

import spkmeans as spka
import spkmeans_module as spkm


def print_2d_list(arr):
    for row in arr:
        for i in range(len(row)):
            row[i] = f'{row[i]:.4f}'
        print(','.join(row))


print_2d_list(spkm.wam_fit([[1.0, 2.0, 3.0],
                            [2.0, 4.0, 6.0],
                            [3.0, 6.0, 8.0]]))


def test_spk():
    print("hello")
