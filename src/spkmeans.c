#include "spkmeans.h"

const char *goal[] = {"wam", "ddg", "lnorm", "jacobi"};

/* Auxiliary methods */
/**
 * @brief free memory of 2d array created by calloc_2d_array()
 *
 * @param arr
 * @param length number of rows
 */
void free_2d(double **arr)
{
    free(arr[0]);
    free(arr);
}
/**
 * @brief creating continious memory block for 2d array
 *
 * @param rows_number number of rows
 * @param col_number number of columns
 * @return double** pointer to arr[0]
 */
double **calloc_2d_array(int rows_number, int cols_number)
{
    int i;
    double *p = calloc(rows_number * cols_number, sizeof(double));

    double **arr = calloc(rows_number, sizeof(double *));
    if (p == NULL || arr == NULL)
        other_error();
    for (i = 0; i < rows_number; i++)
        arr[i] = p + i * cols_number;
    return arr;
}
/**
 * @brief prints "Invalid Input!"
 *
 */
void invalid_input()
{
    printf("Invalid Input!\n");
    exit(0);
}
void other_error()
{
    printf("An Error Has Occurred\n");
    exit(0);
}
/**
 * @brief function to check if the path is valid - ending with .txt or .csv
 *
 * @param s the path to check
 * @return int 1 - valid, 0 - not valid
 */
int is_valid_path(char *s)
{
    char c;
    for (c = *s; c != '\0'; c = *++s)
    {
        if (c == '.')
        {
            s++;
            if (strcmp(s, "txt") == 0 || strcmp(s, "csv") == 0)
                return 1;
            break;
        }
    }
    return 0;
}
/**
 * @brief calculating the norm of the subtraction of two vectors: ||u-v||
 *
 * @param u - array of double
 * @param v - array of double
 * @param dim - dimension of the both arrays
 * @return double - ||u-v||
 */
double subtract_norm(double const *u, double const *v, int dim)
{
    double d = 0;
    int i;
    for (i = 0; i < dim; i++)
        d += (u[i] - v[i]) * (u[i] - v[i]);

    return sqrt(d);
}
/**
 * @brief calculating the dot product of vector u with vector v
 *
 * @param u array of double
 * @param v array of double
 * @param dim
 * @return double
 */
double dot_prod(double *u, double *v, int dim)
{
    double res = 0;
    int i;
    for (i = 0; i < dim; i++)
    {
        res += u[i] * v[i];
    }
    return res;
}
/**
 * @brief Get column j of matrix object as double array
 *
 * @param matrix
 * @param dim
 * @param j
 * @return double*
 */
double *get_column_of_matrix(double **matrix, int dim, int j)
{
    double *col = calloc(dim, sizeof(double));
    int i;
    if (col == NULL)
        other_error();
    for (i = 0; i < dim; i++)
    {
        col[i] = matrix[i][j];
    }

    return col;
}
/**
 * @brief max function between double numbers
 *
 * @param x
 * @param y
 * @return double - max(x,y)
 */
double dmax(double x, double y)
{
    if (x >= y)
        return x;
    return y;
}
/**
 * @brief
 * finding the indexes of the max element in absolute value of the off-diagonal elements of the given matrix.
 *
 * @param matrix
 * @param dim must be >1
 * @return int* array of the indexes of the max element.arr[0] = i, arr[1] = j
 */
int *find_max_element_off_diagonal(double **matrix, int dim)
{
    double max_element;
    int max_i, max_j, i, j;
    int *indexes;
    if (dim == 1)
        return NULL;
    /* max_element is always positive*/
    max_element = fabs(matrix[0][dim - 1]);
    max_i = 0;
    max_j = dim - 1;

    for (i = 0; i < dim; i++)
    {
        for (j = i; j < dim; j++)
        {
            if (i != j)
            {

                if (max_element != dmax(max_element, fabs(matrix[i][j])))
                {
                    max_element = fabs(matrix[i][j]), max_i = i, max_j = j;
                }
            }
        }
    }
    indexes = calloc(2, sizeof(int));
    if (indexes == NULL)
        other_error();
    indexes[0] = max_i, indexes[1] = max_j;
    return indexes;
}
/**
 * @brief retrun the sign of x. if x == 0 then return 1
 *
 * @param x
 * @return double
 */
double sign(double x)
{
    if (x >= 0)
        return 1;
    return -1;
}

/**
 * @brief sum of squares of all off-diagonal elements of given matrix A
 *
 * @param A
 * @param dim
 * @return double
 */
double off(double **A, int dim)
{
    double sum = 0;
    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (i != j)
            {
                sum += pow(A[i][j], 2);
            }
        }
    }
    return sum;
}

/**
 * @brief print square matrix
 *
 * @param arr 2d array
 * @param dim
 */
void print_squared_2d_arr(double **arr, int dim)
{
    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (j == (dim - 1))
                printf("%.4f", arr[i][j]);
            else
                printf("%.4f,", arr[i][j]);
        }
        printf("\n");
    }
}
/**
 * @brief prints 2d array
 *
 * @param arr
 * @param rows_number
 * @param columns_number
 */
void print_2d_arr(double **arr, int rows_number, int columns_number)
{
    int i, j;
    for (i = 0; i < rows_number; i++)
    {
        for (j = 0; j < columns_number; j++)
        {
            if (j == (columns_number - 1))
                printf("%.4f", arr[i][j]);
            else
                printf("%.4f,", arr[i][j]);
        }
        printf("\n");
    }
}

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
void rotation_prod(double **M, int i, int j, double c, double s, int dim)
{
    int temp, k;
    double *i_col, *j_col;
    if (j < i)
    {
        temp = j;
        j = i;
        i = temp;
    }

    i_col = calloc(dim, sizeof(double));

    j_col = calloc(dim, sizeof(double));
    if (j_col == NULL || i_col == NULL)
        other_error();
    i_col[i] = c;
    i_col[j] = -s;
    j_col[i] = s;
    j_col[j] = c;
    for (k = 0; k < dim; k++)
    {
        /* M[k][i] M[k][j]*/
        double k_i = dot_prod(M[k], i_col, dim);
        double k_j = dot_prod(M[k], j_col, dim);
        M[k][i] = k_i, M[k][j] = k_j;
    }
    free(i_col);
    free(j_col);
}
/**
 * @brief summarize vector v (the second) to vector u (the first)
 *
 * @param u first vector
 * @param v second vector to be added to u
 * @param length
 */
void points_sum(double *u, double const *v, int length)
{
    int i;
    for (i = 0; i < length; i++)
        u[i] += v[i];
}

/* Main Methods */
double **wam(double **points, int points_number, int point_dim)
{
    int i, j;

    double **W_matrix = calloc_2d_array(points_number, points_number);
    for (i = 0; i < points_number; i++)
    {
        for (j = 0; j < points_number; j++)
        {
            if (i == j)
                W_matrix[i][j] = 0;
            else
            {
                double calc_norm = subtract_norm(points[i], points[j], point_dim);
                W_matrix[i][j] = exp(-0.5 * calc_norm);
            }
        }
    }

    return W_matrix;
}
double **ddg(double **points, int points_number, int point_dim)
{
    double **W_matrix = wam(points, points_number, point_dim);
    double **D_matrix = calloc_2d_array(points_number, points_number);
    int i, j, k;
    for (i = 0; i < points_number; i++)
    {
        for (j = 0; j < points_number; j++)
        {
            if (i != j)
                D_matrix[i][j] = 0;
            else
            {
                /* summing the i's row of W_matrix */
                for (k = 0; k < points_number; k++)
                {
                    D_matrix[i][j] += W_matrix[i][k];
                }
            }
        }
    }
    free_2d(W_matrix);
    return D_matrix;
}
double **ddg_inverse_square(double **points, int points_number, int point_dim)
{
    double **W_matrix = wam(points, points_number, point_dim);
    double **D_matrix = calloc_2d_array(points_number, points_number);
    int i, j, k;
    for (i = 0; i < points_number; i++)
    {
        for (j = 0; j < points_number; j++)
        {
            if (i != j)
                D_matrix[i][j] = 0;
            else
            {
                /* summing the i's row of W_matrix */
                for (k = 0; k < points_number; k++)
                {
                    D_matrix[i][j] += W_matrix[i][k];
                }
                D_matrix[i][j] = 1 / sqrt(D_matrix[i][j]);
            }
        }
    }
    free_2d(W_matrix);
    return D_matrix;
}
double **lnorm(double **points, int points_number, int point_dim)
{
    double **W_matrix = wam(points, points_number, point_dim);
    double **D_matrix = ddg_inverse_square(points, points_number, point_dim);
    double **lnorm_matrix = calloc_2d_array(points_number, points_number);
    int i, j;
    for (i = 0; i < points_number; i++)
    {
        for (j = 0; j < points_number; j++)
        {
            if (i == j)
            {
                lnorm_matrix[i][i] = 1 - D_matrix[i][i] * D_matrix[i][i] * W_matrix[i][i];
            }
            else
                lnorm_matrix[i][j] = -1 * D_matrix[i][i] * D_matrix[j][j] * W_matrix[i][j];
        }
    }
    free_2d(W_matrix);
    free_2d(D_matrix);
    return lnorm_matrix;
}
double **jacobi(double **A, int dim)
{

    int i, j, k, a, b, r, l;
    double theta, t, c, s;
    double **next_A, **V, **curr_A, **new_V;
    /* initialize V to be unit matrix */
    V = calloc_2d_array(dim, dim);
    for (l = 0; l < dim; l++)
        V[l][l] = 1;

    /* copying A to work with */
    curr_A = calloc_2d_array(dim, dim);
    for (l = 0; l < dim; l++)
        for (k = 0; k < dim; k++)
            curr_A[l][k] = A[l][k];

    for (k = 0; k < 100; k++)
    {
        /* step 1: find the indexes of the max (absolue) element */
        int *indexes = find_max_element_off_diagonal(curr_A, dim);
        i = indexes[0], j = indexes[1];
        free(indexes);
        if (curr_A[i][j] == 0)
        {
            next_A = curr_A;
            break;
        }

        /* step 2: calculate theta -> t -> c -> s */
        theta = (curr_A[j][j] - curr_A[i][i]) / (2 * curr_A[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        /* step 3: calculating V = product of matrix rotations */
        if (k == 0)
        {
            V[i][i] = V[j][j] = c;
            if (j < i)
                V[j][i] = s, V[i][j] = -s;
            else
                V[i][j] = s, V[j][i] = -s;
        }
        else
            rotation_prod(V, i, j, c, s, dim);

        /* setp 4: calculating A': copying A to A' -> updating row i and column j only */
        next_A = calloc_2d_array(dim, dim);

        for (a = 0; a < dim; a++)
            for (b = 0; b < dim; b++)
                next_A[a][b] = curr_A[a][b];

        for (r = 0; r < dim; r++)
        {
            if (r != i && r != j)
            {
                next_A[r][i] = next_A[i][r] = c * curr_A[r][i] - s * curr_A[r][j];

                next_A[r][j] = next_A[j][r] = c * curr_A[r][j] + s * curr_A[r][i];
            }
        }
        next_A[i][i] = pow(c, 2) * curr_A[i][i] + pow(s, 2) * curr_A[j][j] - 2 * s * c * curr_A[i][j];
        next_A[j][j] = pow(s, 2) * curr_A[i][i] + pow(c, 2) * curr_A[j][j] + 2 * s * c * curr_A[i][j];
        next_A[i][j] = next_A[j][i] = 0;

        /* step 5: check if we reach the required convergence - 1*10^-5 */
        if (off(curr_A, dim) - off(next_A, dim) <= 0.00001)
        {
            free_2d(curr_A);
            break;
        }
        /* step 6: free memory of A and updating A to be next_A */
        free_2d(curr_A);
        curr_A = next_A;
    }

    /* step 7: return new matrix where its first row are the eiganvalues and the rest is matrix V */
    new_V = calloc_2d_array(dim + 1, dim);
    for (i = 0; i < dim + 1; i++)
    {
        for (j = 0; j < dim && i > 0; j++)
        {
            new_V[0][j] = next_A[j][j];
            if (V[i - 1][j] == 0)
            {
                V[i - 1][j] = 0;
            }

            new_V[i][j] = V[i - 1][j];
        }
    }
    free_2d(V);
    free_2d(next_A);
    return new_V;
}
double **k_means(int const K, int const max_iter, double **points, double **centroids, int const points_length, int const point_length, double eps)
{
    int k;
    int i;
    int j;
    int l;
    int all_diffs_are_small;
    int *points_to_centroids;
    int *diff_small;
    int new_centroids_index;
    double *avg_centroid;
    double counter;

    points_to_centroids = calloc(points_length, sizeof(int));
    diff_small = calloc(K, sizeof(int));
    if (diff_small == NULL || points_to_centroids == NULL)
        other_error();

    for (i = 0; i < K; i++)
    {
        diff_small[i] = 1;
    }

    for (l = 0; l < max_iter; l++)
    {
        for (i = 0; i < points_length; i++)
        {
            new_centroids_index = points_to_centroids[i];
            for (j = 0; j < K; j++)
            {
                if (subtract_norm(points[i], centroids[j], point_length) < subtract_norm(points[i], centroids[new_centroids_index], point_length))
                {
                    new_centroids_index = j;
                }
            }
            points_to_centroids[i] = new_centroids_index;
        }

        for (i = 0; i < K; i++)
        {
            diff_small[i] = 0;
            avg_centroid = calloc(point_length, sizeof(double));
            if (avg_centroid == NULL)
                other_error();
            counter = 0;
            for (j = 0; j < points_length; j++)
            {
                if (points_to_centroids[j] == i)
                {
                    points_sum(avg_centroid, points[j], point_length);
                    counter++;
                }
            }

            if (counter > 0)
            {
                for (k = 0; k < point_length; k++)
                    avg_centroid[k] = avg_centroid[k] / counter;

                if (subtract_norm(centroids[i], avg_centroid, point_length) < eps)
                {
                    diff_small[i] = 1;
                }
                free(centroids[i]);
                centroids[i] = avg_centroid;
            }
        }
        all_diffs_are_small = 0;
        for (k = 0; k < K; k++)
        {
            all_diffs_are_small += diff_small[k];
        }
        if (all_diffs_are_small == K)
        {
            break;
        }
    }
    free(points_to_centroids);
    free(diff_small);
    return centroids;
}

int main(int argc, char **argv)
{
    char *input_goal;
    char *input_path;
    int i, j;
    int point_dim = 0, points_number = 0;
    char c;
    double **points;
    FILE *fptr;
    if (argc != 3)
        invalid_input();
    input_goal = argv[1];
    input_path = argv[2];
    if (!is_valid_path(input_path))
        invalid_input();

    /* reading input from file */
    fptr = fopen(input_path, "r");
    if (fptr == NULL)
    {
        other_error();
    }

    /* counting point size */
    while ((c = fgetc(fptr)) != EOF)
    {
        double point;
        if (c == '\n')
        {
            rewind(fptr);
            break;
        }
        point_dim++;
        fscanf(fptr, "%lf", &point);
    }

    /* counting number of points */
    while ((c = fgetc(fptr)) != EOF)
        if (c == '\n')
            points_number++;

    rewind(fptr);

    points = calloc_2d_array(points_number, point_dim);

    /*  reading the points from the file into array */
    for (i = 0; i < points_number; i++)
    {
        /* points[i] = p + i * point_dim; */
        for (j = 0; j < point_dim; j++)
        {
            fscanf(fptr, "%lf", &points[i][j]);
            fgetc(fptr);
        }
    }
    fclose(fptr);

    if (!strcmp(input_goal, goal[0]))
    {
        double **matrix = wam(points, points_number, point_dim);
        print_squared_2d_arr(matrix, points_number);
        free_2d(matrix);
    }
    else if (!strcmp(input_goal, goal[1]))
    {
        double **matrix = ddg(points, points_number, point_dim);
        print_squared_2d_arr(matrix, points_number);
        free_2d(matrix);
    }
    else if (!strcmp(input_goal, goal[2]))
    {
        double **matrix = lnorm(points, points_number, point_dim);
        print_squared_2d_arr(matrix, points_number);
        free_2d(matrix);
    }
    else if (!strcmp(input_goal, goal[3]))
    {
        /* we assume the input is a symmetric matrix */
        double **matrix = jacobi(points, points_number);

        print_2d_arr(matrix, points_number + 1, points_number);
        free_2d(matrix);
    }
    else
        invalid_input();

    free_2d(points);

    return 0;
}
