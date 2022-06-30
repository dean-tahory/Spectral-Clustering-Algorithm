import pytest
import spkmeans

import pytest
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

# https://www.geeksforgeeks.org/ml-spectral-clustering/
# https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html#sklearn.cluster.SpectralClustering.fit
# https://docs.pytest.org/en/7.1.x/getting-started.html#getstarted


# Loading the data
points_table = pd.read_csv('input.txt', header=None)
points = points_table.to_numpy()
print(points_table.head())

sc = SpectralClustering(n_clusters=2, affinity='precomputed')
w_matrix = spkmeans.wam_fit(points.tolist(), points.shape[0])
w_matrix_table = pd.DataFrame(w_matrix)
print(w_matrix_table.head())

# sc.fit_predict(w_matrix)
centroids = sc.fit_predict(points_table)
print(centroids)
