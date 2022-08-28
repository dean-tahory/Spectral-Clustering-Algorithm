
double **calloc_2d_array(int rows_number, int cols_number);
void other_error();
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
 * @brief return centroids acording to spectral clustering algorithm as array of double arrays which represent the points.
 *
 * @param K number of centroids
 * @param points 2d list of points to cluster
 * @param centroids initial centroids
 * @return list of centroids points (2d array)
 */
double **spk(int const K, double **points, double **centroids);

/**
 * @brief
 *
 * @param A
 * @param dim dimension of matrix A
 * @return double** 2d array where the first row is the eiganvalues and the rest is matrix V which its columns are eiganvectors of A.
 *         more precisely, the column below eignvalue at the first row is its eiganvector
 */
double **jacobi(double **A, int dim);

double **k_means(int const K, int const max_iter, double **points, double **centroids, int const points_length, int const point_length, double eps);

/**
 * @brief free memory of 2d array created by calloc_2d_array()
 *
 * @param arr
 * @param length number of rows
 */
void free_2d(double **arr);

/**
 * @brief creating continious memory block for 2d array
 *
 * @param rows_number number of rows
 * @param col_number number of columns
 * @return double** pointer to arr[0]
 */
double **calloc_2d_array(int rows_number, int cols_number);

double norm(double const *u, double const *v, int dim);
double dot_prod(double *u, double *v, int dim);
double *get_column_of_matrix(double **matrix, int dim, int j);
double dmax(double x, double y);
/**
 * @brief
 * finding the indexes of the max element in absolute value of the off-diagonal elements of the given matrix.
 *
 * @param matrix
 * @param dim must be >1
 * @return int* array of the indexes of the max element.arr[0] = i, arr[1] = j
 */
int *find_max_element_off_diagonal(double **matrix, int dim);
/**
 * @brief retrun the sign of x. if x == 0 then return 1
 *
 * @param x
 * @return double
 */
double sign(double x);

/**
 * @brief sum of squares of all off-diagonal elements of given matrix A
 *
 * @param A
 * @param dim
 * @return double
 */
double off(double **A, int dim);

/**
 * @brief print square matrix
 *
 * @param arr 2d array
 * @param dim
 */
void print_squared_2d_arr(double **arr, int dim);

void print_2d_arr(double **arr, int rows_number, int columns_number);

/**
 * @brief efficient implementation of matrix product where the second matrix is a rotation matrix
 * denote P =  rotation matrix - unit matrix where at P[i][i] = c, P[i][j] = s, P[j][j] = -s, P[j][j] = c
 * @param M the matrix to be rotated
 * @param i
 * @param j
 * @param c value of P[i][i] and P[j][j]
 * @param s value of P[i][j] and -P[j][i]
 * @param dim
 */
void rotation_prod(double **M, int i, int j, double c, double s, int dim);
/* sum the second vector v into u */
void points_sum(double *u, double const *v, int length);
