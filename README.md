# Spectral Clustering Algorithm

## Overview

As part of my Software Project course at Tel Aviv University, I've been asked to implement the normalized version of spectral clustering algorithm using c and export it as a python extension module.

## Normalized Spectral Clustering

This is the Normalized Spectral Clustering algorithm based on [[1,2]](https://github.com/dean-tahory/Spectral-Clustering-Algorithm/new/main?readme=1#references). Given a set of $n$ points $X = x_1,x_2,...x_N\in R^d$ the algorithm is:

1. Form the weighted adjacency matrix $W$ from $X$
2. Compute the normalized graph Laplacian $L_{norm}$
3. If given $k$ is equal to $0$ we determine it using the [eigangape heuristic]()
   and obtain the largest $k$ eiganvalues $u_1,...u_k$ of $L_{norm}$
4. let $U\in R^{n\times k}$ be the matrix containing the vectors $u_1,...u_k$ as columns
5. Form the matrix $T\in R^{n\times k}$ from $U$ by renormalizing each of $U$'s rows to have unit length, that is set $t_{ij}= \frac{u_{ij}}{(\sum_j{u_{ij}^2})^{\frac{1}{2}}}$
6. Treating each row of $T$ as point in $R^k$, cluster them into $k$ clusters via the K-means algorithm

### The Eigengap Heuristic

In order to determine the number of clusters $k$, we will use eigengap heuristic as follow:
let $(δ_i)_{i=1,...,n−1}$ be the eigengap for the decreasing ordered eigenvalues $λ1 ≥ ... ≥ λn ≥ 0$ of $L_{norm}$, defined as:
$$δ_i = |λ_i − λ_{i+1}|$$
The eigengap measure could indicate the number of clusters $k$ through the following estimation:
$$k = argmax_i(δ_i), i = 1,...,\lfloor{\frac{n}{2}}\rfloor$$
In case of equality in the argmax of some eigengaps, we use the lowest index.

## Features

- implemented Jacobi eiganvalue algorithm which calculates the eiganvalues and eiganvectors of real symmetric matrix (diagonalization)
- implemented k-means clustering algorithm
- used `pytest` for unit testing

## How to install and run the project

### Installing

1. Clone the project repository, then on command line run: `cd project`
2. Install the dependencies: `pipenv install`
3. activate the Pipenv shell: `pipenv shell`
4. create the python extension module

```
python3 setup.py build_ext --inplace
```

5. complile the c files using the `comp.sh` shell script:

```
bash comp.sh
```

now you ready to run the algorithm (`spkmeans.py`).

### Running

the format of the CMD arguments is as follows:

1. _k_ (int, $< N$) : Number of required clusters. If equal 0, we use the eigangap heuristic. MUST be passed for all goals.
2. _goal_ (enum): Can get the following values:
   - _spk_: Perform full spectral kmeans
   - _wam_: Calculate and output the Weighted Adjacency Matrix
   - _ddg_: Calculate and output the Diagonal Degree Matrix
   - _lnorm_: Calculate and output the Normalized Graph Laplacian
   - _jacobi_: Calculate and output the eigenvalues and eigenvectors
3. _file_name_ (`.txt` or `.csv`): The path to the Input file, it will contain N data points for all above goals except Jacobi, in case the goal is Jacobi the input is a symmetric matrix, the file extension is .txt or .csv.
   - the .txt file must have new line at the end.
   - the points coordiantes are separated by comma ','
   - each point is written in its own line

Example:

```
python3 spkmeans.py 3 spk input_1.txt
```

### Output format

- Case of ’_spk_’: The first line will be the indices of the observations chosen by the K-means++ algorithm as the initial centroids. We refer to the first observation index as 0, the second as 1 and so on, up until N - 1. The second line onward will be the calculated final centroids from the K-means algorithm, separated by a comma, such that each centroid is in a line of its own.
- Case of ’_Jacobi_’: The first line will be the eigenvalues, second line onward will be the corresponding eigenvectors (printed as columns).
- Else: output the required matrix separated by a comma, such that each row is in a line of its own.

## References

- [1] Andrew Ng, Michael Jordan, and Yair Weiss. On spectral clustering: Analysis and an algorithm. Advances in neural information processing systems, 14:849–856, 2001.
- [2] Ulrike Von Luxburg. A tutorial on spectral clustering. Statistics and computing, 17(4):395–416, 2007.
