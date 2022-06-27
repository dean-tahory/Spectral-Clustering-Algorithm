#include "spkmeans.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// TODO remove 'test label
const char *goal[] = {"wam", "ddg", "lnorm", "jacobi", "test"};

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
    int res = 0;
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
            C[i][j] = dot_prod(A[i], get_column_of_matrix(B, dim, j), dim);
        }
    }
    return C;
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
    return lnorm_matrix;
}

// testing methods
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
    // TODO it tested but it's good to add here test
}
void test_matrix_prod()
{
    double my_A[3][3] = {
        {2, 1, 2},
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
    double my_B[3][3] = {
        {2, 1, 2},
        {5, 9, 4},
        {7, 5, 8}};
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
void test()
{
    test_matrix_prod();
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