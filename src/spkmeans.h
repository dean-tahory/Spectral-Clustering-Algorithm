#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief allocating continuous memory block for 2d array
 *
 * @param rows_number
 * @param cols_number
 * @return double**
 */
double **calloc_2d_array(int rows_number, int cols_number);
/**
 * @brief prints "An Error Has Occurred"
 *
 */
void other_error(void);
/**
 * @brief free the continuous memmory block which allocated with calloc_2d_array()
 *
 * @param arr
 */
void free_2d(double **arr);
/**
 * @brief calculating the weighted adjacency matrix
 *
 * @param points array of points represented as array
 * @param points_number number of points
 * @param point_dim dimension of point
 * @return 2d array
 */
double **wam(double **points, int points_number, int point_dim);

/**
 * @brief calculating the diagonal degree matrix
 *
 * @param points array of points represented as array
 * @param points_number number of points
 * @param point_dim dimension of point
 * @return 2d array
 */
double **ddg(double **points, int points_number, int point_dim);

/**
 * @brief calculating the normalized graph laplacian
 *
 * @param points array of points represented as array
 * @param points_number number of points
 * @param point_dim dimension of point
 * @return 2d array
 */
double **lnorm(double **points, int points_number, int point_dim);

/**
 * @brief
 *
 * @param A matrix
 * @param dim dimension of matrix A
 * @return double** 2d array where the first row is the eiganvalues and the rest is matrix V which its columns are eiganvectors of A.
 *         more precisely, the column below eignvalue at the first row is its eiganvector
 */
double **jacobi(double **A, int dim);
/**
 * @brief K-means algorithm which clusters the given points to K centoids and return them
 *
 * @param K number of desired clusters
 * @param max_iter maximum of iterations to allow
 * @param points points to cluster
 * @param centroids initial centroids for the algorithm
 * @param points_length number of points
 * @param point_length dimension of point
 * @param eps epsilon - for desired precision
 * @return double** array of centroids
 */
double **k_means(int const K, int const max_iter, double **points, double **centroids, int const points_length, int const point_length, double eps);
