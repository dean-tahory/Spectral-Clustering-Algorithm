

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
 * @brief return centroids acording to spectral clustering algorithm as array of double arrays which represent the points.
 *
 * @param K number of centroids
 * @param points 2d list of points to cluster
 * @param centroids initial centroids
 * @return list of centroids points (2d array)
 */
double **spk(int const K, double **points, double **centroids);