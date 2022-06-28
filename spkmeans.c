#include "spkmeans.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// TODO remove 'test label
const char *goal[] = {"wam", "ddg", "lnorm", "jacobi", "test"};
// TODO move this method to test section
void print_arr(double *arr, int length)
{
    for (int i = 0; i < length; i++)
    {
        if (i == (length - 1))
            printf("%.4lf", arr[i]);
        else
            printf("%.4lf,", arr[i]);
    }
}

/* Auxiliary methods */
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
int is_valid_path(char *s)
{
    int dot_flag = 0;
    for (char c = *s; c != '\0'; c = *++s)
    {
        if (c == '.')
        {
            dot_flag = 1;
            if (strcmp(++s, "txt") != 0)
                return 0;
        }
    }
    if (!dot_flag)
        return 0;
    return 1;
}
double norm(double const *u, double const *v, int dim)
{
    double d = 0;
    int i;
    for (i = 0; i < dim; i++)
        d += (u[i] - v[i]) * (u[i] - v[i]);

    return sqrt(d);
}
double dot_prod(double *u, double *v, int dim)
{
    double res = 0;
    for (int i = 0; i < dim; i++)
    {
        res += u[i] * v[i];
    }
    return res;
}
double *get_column_of_matrix(double **matrix, int dim, int j)
{
    double *col = calloc(dim, sizeof(double));
    for (int i = 0; i < dim; i++)
    {
        col[i] = matrix[i][j];
    }

    return col;
}
/**
 * @brief multiply square matrices
 *
 * @param A 2d array
 * @param B 2d array
 * @param dim  dim of the matrices
 */
// TODO remove this method if it's not used
double **matrix_prod(double **A, double **B, int dim)
{
    double **C = calloc(dim, sizeof(double *));
    for (int i = 0; i < dim; i++)
    {
        C[i] = calloc(dim, sizeof(double));
        for (int j = 0; j < dim; j++)
        {
            // printf("A[%d]: ", i);
            // print_arr(A[i], dim);
            // printf("\nB[%d]: ", j);
            // print_arr(get_column_of_matrix(B, dim, j), dim);
            // printf("\n");
            C[i][j] = dot_prod(A[i], get_column_of_matrix(B, dim, j), dim);
            // printf("c[%d][%d]: %.4lf\n", i, j, C[i][j]);
        }
    }
    return C;
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
    if (dim == 1)
        return NULL;
    // max_element is always positive
    double max_element = fabs(matrix[0][dim - 1]);
    int *max_index = calloc(2, sizeof(int));
    int max_i = 0, max_j = dim - 1;
    for (int i = 0; i < dim; i++)
    {
        for (int j = i; j < dim; j++)
        {
            if (i != j)
            {
                if (max_element != fmax(max_element, fabs(matrix[i][j])))
                {
                    max_element = fabs(matrix[i][j]), max_i = i, max_j = j;
                }
            }
        }
    }
    int *indexes = calloc(2, sizeof(int));
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
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
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
void print_2d_arr(double **arr, int dim)
{
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            if (j == (dim - 1))
                printf("%.4lf", arr[i][j]);
            else
                printf("%.4lf,", arr[i][j]);
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
    if (j < i) // swap, we want p,p to be the first element in the diagonal where it equals to a
    {
        int temp = j;
        j = i;
        i = temp;
    }

    double *i_col = calloc(dim, sizeof(double));

    double *j_col = calloc(dim, sizeof(double));

    i_col[i] = c;
    i_col[j] = -s;
    j_col[i] = s;
    j_col[j] = c;

    for (int k = 0; k < dim; k++)
    {
        // M[k][i] M[k][j]
        double k_i = dot_prod(M[k], i_col, dim);
        double k_j = dot_prod(M[k], j_col, dim);
        M[k][i] = k_i, M[k][j] = k_j;
    }
}

/* Main Methods */
double **wam(double **points, int points_number, int point_dim)
{
    double **W_matrix = calloc(points_number, sizeof(double *));
    for (int i = 0; i < points_number; i++)
    {
        W_matrix[i] = calloc(points_number, sizeof(double));
        for (int j = 0; j < points_number; j++)
        {
            if (i == j)
                W_matrix[i][j] = 0;
            else
                W_matrix[i][j] = exp(-0.5 * norm(points[i], points[j], point_dim));
        }
    }

    return W_matrix;
}
double **ddg(double **points, int points_number, int point_dim)
{
    double **W_matrix = wam(points, points_number, point_dim);
    double **D_matrix = calloc(points_number, sizeof(double *));
    for (int i = 0; i < points_number; i++)
    {
        D_matrix[i] = calloc(points_number, sizeof(double));
        for (int j = 0; j < points_number; j++)
        {
            if (i != j)
                D_matrix[i][j] = 0;
            else
            {
                // summing the i's row of W_matrix
                for (int k = 0; k < points_number; k++)
                {
                    D_matrix[i][j] += W_matrix[i][k];
                }
                D_matrix[i][j] = 1 / sqrt(D_matrix[i][j]);
            }
        }
    }
    free(W_matrix);
    return D_matrix;
}
double **lnorm(double **points, int points_number, int point_dim)
{
    double **W_matrix = wam(points, points_number, point_dim);
    double **D_matrix = ddg(points, points_number, point_dim);
    double **lnorm_matrix = calloc(points_number, sizeof(double *));
    for (int i = 0; i < points_number; i++)
    {
        lnorm_matrix[i] = calloc(points_number, sizeof(double));
        for (int j = 0; j < points_number; j++)
        {
            if (i == j)
            {
                lnorm_matrix[i][i] = 1 - D_matrix[i][i] * D_matrix[i][i] * W_matrix[i][i];
            }
            else
                lnorm_matrix[i][j] = D_matrix[i][i] * D_matrix[j][j] * W_matrix[i][j];
        }
    }
    free(W_matrix);
    free(D_matrix);
    return lnorm_matrix;
}
double **jacobi(double **A, int dim)
{

    int i, j;
    double theta, t, c, s;
    double **next_A;
    // initialize V to be unit matrix
    double **V = calloc(dim, sizeof(double *));
    for (int l = 0; l < dim; l++)
    {
        V[l] = calloc(dim, sizeof(double));
        V[l][l] = 1;
    }

    for (int k = 0; k < 100; k++)
    {
        // step 1: find the indexes of the max (absolue) element
        int *indexes = find_max_element_off_diagonal(A, dim);
        i = indexes[0], j = indexes[1];

        // step 2: calculate theta -> t -> c -> s
        theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
        t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
        c = 1 / (sqrt(pow(t, 2) + 1));
        s = t * c;

        // step 3: calculating V = product of matrix rotations
        if (k == 0)
        {
            V[i][i] = V[j][j] = c;
            if (j < i)
                V[j][i] = s, V[i][j] = -s;
            else
                V[i][j] = s, V[j][i] = -s;
        }
        else
        {
            /* using matrix product */
            // double **P = calloc(dim, sizeof(double *));
            // for (int l = 0; l < dim; l++)
            // {
            //     P[l] = calloc(dim, sizeof(double));
            //     P[l][l] = 1;
            // }
            // P[i][i] = P[j][j] = c;
            // if (j < i)
            //     P[j][i] = s, P[i][j] = -s;
            // else
            //     P[i][j] = s, P[j][i] = -s;
            // V = matrix_prod(V, P, dim);

            rotation_prod(V, i, j, c, s, dim);
        }

        // setp 4: calculating A': copying A to A' -> updating row i and column j only
        next_A = calloc(dim, sizeof(double *));
        for (int a = 0; a < dim; a++)
        {
            next_A[a] = calloc(dim, sizeof(double));
            for (int b = 0; b < dim; b++)
            {
                next_A[a][b] = A[a][b];
            }
        }

        for (int r = 0; r < dim; r++)
        {
            if (r != i && r != j)
            {
                next_A[r][i] = next_A[i][r] = c * A[r][i] - s * A[r][j];

                next_A[r][j] = next_A[j][r] = c * A[r][j] + s * A[r][i];
            }
        }
        next_A[i][i] = pow(c, 2) * A[i][i] + pow(s, 2) * A[j][j] - 2 * s * c * A[i][j];
        next_A[j][j] = pow(s, 2) * A[i][i] + pow(c, 2) * A[j][j] + 2 * s * c * A[i][j];
        next_A[i][j] = next_A[j][i] = 0;

        // step 5: check if we reach the required convergence - 1*10^-5
        if (off(A, dim) - off(next_A, dim) <= pow(10, -5))
        {
            break;
        }
        // step 6: free memory of A and updating A to be next_A
        free(A);
        A = next_A;
    }

    // step 7: print the eiganvalues
    for (int l = 0; l < dim; l++)
    {
        if (l == dim - 1)
            printf("%.4lf\n", next_A[l][l]);
        else
            printf("%.4lf,", next_A[l][l]);
    }

    return V;
}

// testing methods
void test_get_column_of_matrix()
{
    double my_mat[3][3] = {
        {2, 1, 2},
        {3, 3, 4},
        {2, 5, 8}};
    double **mat = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        mat[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            mat[i][j] = my_mat[i][j];
        }
    }
    print_arr(get_column_of_matrix(mat, 3, 2), 3);
    printf("\n");
}
void test_dot_prod()
{
}
void test_matrix_prod()
{
    double my_A[3][3] = {
        {1, 0, 0},
        {0, 0.8112, 0.5847},
        {0, -0.5847, 0.8112}};
    double **A = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        A[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = my_A[i][j];
        }
    }
    double my_B[3][3] = {
        {1, 0, 0},
        {0, 0.8112, 0.5847},
        {0, -0.5847, 0.8112}};
    double **B = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        B[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            B[i][j] = my_B[i][j];
        }
    }
    print_2d_arr(matrix_prod(A, B, 3), 3);
}
void test_find_max_element_off_diagonal()
{
    double my_A[3][3] = {
        {2, 1, -50},
        {3, 3, 4},
        {2, 5, 8}};
    double **A = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        A[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = my_A[i][j];
        }
    }
    int *index = find_max_element_off_diagonal(A, 3);
    printf("i:%d,j:%d", index[0], index[1]);
}
void test_off()
{
    double my_A[3][3] = {
        {1, 2, -3},
        {2, 4, 6},
        {3, 6, 8}};
    double **A = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        A[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = my_A[i][j];
        }
    }
    printf("%.4lf\n", off(A, 3));
}
void test_jacobi()
{
    // usefull link to test - https://tiffzhang.com/jacobi/ + http://www.math.utoledo.edu/~codenth/Linear_Algebra/Calculators/jacobi_algorithm.html
    double my_A[3][3] = {
        {1, 2, 3},
        {2, 4, 6},
        {3, 6, 8}};
    double **A = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        A[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = my_A[i][j];
        }
    }
    print_2d_arr(jacobi(A, 3), 3);
}
void test_rotation_prod()
{
    // TODO test with higher dimension matrix
    // double my_A[4][4] = {
    //     {1, 2, 3, 4},
    //     {5, 6, 7, 8},
    //     {9, 10, 11, 12},
    //     {13, 14, 15, 16}};
    // double **A = calloc(4, sizeof(double *));
    // for (int i = 0; i < 4; i++)
    // {
    //     A[i] = calloc(4, sizeof(double));
    //     for (int j = 0; j < 4; j++)
    //     {
    //         A[i][j] = my_A[i][j];
    //     }
    // }
    // double my_B[4][4] = {
    //     {1, 0, 0, 0},
    //     {0, 2, 0, 4},
    //     {0, 0, 1, 0},
    //     {0, -4, 0, 2}};
    // double **B = calloc(4, sizeof(double *));
    // for (int i = 0; i < 4; i++)
    // {
    //     B[i] = calloc(4, sizeof(double));
    //     for (int j = 0; j < 4; j++)
    //     {
    //         B[i][j] = my_B[i][j];
    //     }
    // }
    double my_A[3][3] = {
        {1.0000, 0.0000, 0.0000},
        {0.0000, 0.8112, 0.5847},
        {0.0000, -0.5847, 0.8112}};
    double **A = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        A[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = my_A[i][j];
        }
    }
    double my_B[3][3] = {
        {1.0000, 0.0000, 0.0000},
        {0.0000, 0.8112, 0.5847},
        {0.0000, -0.5847, 0.8112}};
    double **B = calloc(3, sizeof(double *));
    for (int i = 0; i < 3; i++)
    {
        B[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; j++)
        {
            B[i][j] = my_B[i][j];
        }
    }
    rotation_prod(A, 1, 2, 0.8112, 0.5847, 3);
    print_2d_arr(A, 3);
}

void test()
{
    test_jacobi();
    // test_rotation_prod();
}
int main(int argc, char **argv)
{
    char *input_goal;
    char *input_path;

    if (argc != 3)
        invalid_input();
    input_goal = argv[1];
    input_path = argv[2];
    if (!is_valid_path(input_path))
        invalid_input();

    /* reading input from file */
    FILE *fptr = fopen(input_path, "r");
    if (fptr == NULL)
    {
        other_error();
    }

    int point_dim = 0;
    int points_number = 0;
    char c;

    // counting point size
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

    // counting number of points
    while ((c = fgetc(fptr)) != EOF)
        if (c == '\n')
            points_number++;

    rewind(fptr);

    double **points = calloc(points_number, sizeof(double *));
    // reading the points from the file into array
    for (int i = 0; i < points_number; i++)
    {
        points[i] = calloc(point_dim, sizeof(double));
        for (int j = 0; j < point_dim; j++)
        {
            fscanf(fptr, "%lf", &points[i][j]);
            fgetc(fptr);
        }
    }
    fclose(fptr);

    if (!strcmp(input_goal, goal[0]))
    {
        print_2d_arr(wam(points, points_number, point_dim), points_number);
    }
    if (!strcmp(input_goal, goal[1]))
    {
        print_2d_arr(ddg(points, points_number, point_dim), points_number);
    }

    // TODO remove 'test' label case
    if (!strcmp(input_goal, goal[4]))
    {
        test();
    }

    else
        invalid_input();

    // free memory
    // for (int i = 0; i < K; i++)
    //     free(centroids[i]);
    // free(centroids);
    for (int i = 0; i < points_number; i++)
        free(points[i]);
    free(points);

    return 0;
}