#include "spkmeans.c"
#include <Python.h>
#include <assert.h>
// method to convert python list of lists to 2d array in C
double **python_matrix_to_c(PyObject *_list, int points_number, int point_dim)
{
    PyObject *point, *point_value;
    Py_ssize_t i, j;

    double **points = calloc_2d_array(points_number, point_dim);
    for (i = 0; i < points_number; i++)
    {
        point = PyList_GetItem(_list, i);
        if (!PyList_Check(point))
        { /* We only adding lists of points */
            continue;
        }

        for (j = 0; j < point_dim; j++)
        {
            point_value = PyList_GetItem(point, j);
            if (!PyFloat_Check(point_value))
                continue; /* Skip non-float values */

            points[i][j] = PyFloat_AS_DOUBLE(point_value); /* Convert a Python float object to double */
            if (points[i][j] == -1 && PyErr_Occurred())
            {
                /* Float too big to fit in a C double, bail out */
                other_error();
            }
        }
    }
    return points;
}

double **matrix_to_not_continuous_matrix(double **matrix, int m, int n)
{
    double **arr = calloc(m, sizeof(double *));
    for (int i = 0; i < m; i++)
    {
        arr[i] = calloc(n, sizeof(double));
        for (int j = 0; j < n; j++)
        {
            *(arr[i] + j) = matrix[i][j];
        }
    }
    free_2d(matrix);
    return arr;
}

//  our main module API method - getting args from python object, converting them to C objects and work on them.
static PyObject *kmeans_fit(PyObject *self, PyObject *args)
{
    int K, max_iter, points_number, point_dim, i;
    double eps;
    PyObject *points_list, *centroids_list;

    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "iiOOiid:fit", &K, &max_iter, &points_list, &centroids_list, &points_number, &point_dim, &eps))
    {
        return NULL;
    }
    /* check if pointslits is of type list */
    if (!PyList_Check(points_list))
        return NULL;
    double **points = python_matrix_to_c(points_list, points_number, point_dim);

    // the spk function requires array of arrays where each array is allocated on its on
    double **initial_centroids = matrix_to_not_continuous_matrix(python_matrix_to_c(centroids_list, K, point_dim), K, point_dim);

    // calling the kmeans C algorithm (new version that also gets initial centroids from the python module)
    double **centroids = k_means(K, max_iter, points, initial_centroids, points_number, point_dim, eps);

    // converting C 2d array to Python list of lists.
    PyObject *new_centroids = PyList_New(0);
    if (new_centroids == NULL)
        other_error();

    for (i = 0; i < K; i++)
    {
        PyObject *new_centroid_point = PyList_New(0);
        if (new_centroid_point == NULL)
            other_error();
        Py_ssize_t j;
        for (j = 0; j < point_dim; j++)
        {
            // not printing centroifs[i][j] because its not continiuous
            PyObject *float_point = PyFloat_FromDouble(*(centroids[i] + j));
            PyList_Append(new_centroid_point, float_point);
        }
        if (PyList_Append(new_centroids, new_centroid_point) == -1)
            other_error();
    }

    free_2d(points);
    free_2d(centroids);

    return new_centroids;
}

static PyObject *wam_fit(PyObject *self, PyObject *args)
{
    int points_number, point_dim, i, j;
    PyObject *points_list;
    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "O:wam_fit", &points_list))
    {
        return NULL;
    }

    /* check if pointslits is of type list */
    if (!PyList_Check(points_list))
        return NULL;
    /* Get the size of it and build the output list */
    points_number = PyList_Size(points_list); /*  Same as in Python len(_list)  */
    point_dim = PyList_Size(PyList_GetItem(points_list, 0));

    double **points = python_matrix_to_c(points_list, points_number, point_dim);
    double **matrix = wam(points, points_number, point_dim);

    // converting C 2d array to Python list of lists.
    PyObject *py_matrix = PyList_New(0);
    if (py_matrix == NULL)
        other_error();

    for (i = 0; i < points_number; i++)
    {
        PyObject *row_i = PyList_New(0);
        if (row_i == NULL)
            other_error();
        for (j = 0; j < points_number; j++)
        {
            PyObject *float_point = PyFloat_FromDouble((double)matrix[i][j]);
            PyList_Append(row_i, float_point);
        }
        if (PyList_Append(py_matrix, row_i) == -1)
            other_error();
    }

    // free memory

    free_2d(matrix);
    free_2d(points);
    return py_matrix;
}

static PyObject *ddg_fit(PyObject *self, PyObject *args)
{
    int points_number, point_dim, i, j;
    PyObject *points_list;
    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "O:wam_fit", &points_list))
    {
        return NULL;
    }

    /* check if pointslits is of type list */
    if (!PyList_Check(points_list))
        return NULL;
    /* Get the size of it and build the output list */
    points_number = PyList_Size(points_list); /*  Same as in Python len(_list)  */
    point_dim = PyList_Size(PyList_GetItem(points_list, 0));

    double **points = python_matrix_to_c(points_list, points_number, point_dim);
    double **matrix = ddg(points, points_number, point_dim);

    // converting C 2d array to Python list of lists.
    PyObject *py_matrix = PyList_New(0);
    if (py_matrix == NULL)
        other_error();

    for (i = 0; i < points_number; i++)
    {
        PyObject *row_i = PyList_New(0);
        if (row_i == NULL)
            other_error();
        for (j = 0; j < points_number; j++)
        {
            PyObject *float_point = PyFloat_FromDouble((double)matrix[i][j]);
            PyList_Append(row_i, float_point);
        }
        if (PyList_Append(py_matrix, row_i) == -1)
            other_error();
    }

    // free memory

    free_2d(matrix);
    free_2d(points);
    return py_matrix;
}

static PyObject *lnorm_fit(PyObject *self, PyObject *args)
{
    int points_number, point_dim, i, j;
    PyObject *points_list;
    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "O:wam_fit", &points_list))
    {
        return NULL;
    }

    /* check if pointslits is of type list */
    if (!PyList_Check(points_list))
        return NULL;
    /* Get the size of it and build the output list */
    points_number = PyList_Size(points_list); /*  Same as in Python len(_list)  */
    point_dim = PyList_Size(PyList_GetItem(points_list, 0));

    double **points = python_matrix_to_c(points_list, points_number, point_dim);
    double **matrix = lnorm(points, points_number, point_dim);

    // converting C 2d array to Python list of lists.
    PyObject *py_matrix = PyList_New(0);
    if (py_matrix == NULL)
        other_error();

    for (i = 0; i < points_number; i++)
    {
        PyObject *row_i = PyList_New(0);
        if (row_i == NULL)
            other_error();
        for (j = 0; j < points_number; j++)
        {
            PyObject *float_point = PyFloat_FromDouble((double)matrix[i][j]);
            PyList_Append(row_i, float_point);
        }
        if (PyList_Append(py_matrix, row_i) == -1)
            other_error();
    }

    // free memory

    free_2d(matrix);
    free_2d(points);
    return py_matrix;
}

static PyObject *jacobi_fit(PyObject *self, PyObject *args)
{
    int points_number, point_dim, i, j;
    PyObject *points_list;
    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "O:wam_fit", &points_list))
    {
        return NULL;
    }

    /* check if pointslits is of type list */
    if (!PyList_Check(points_list))
        return NULL;
    /* Get the size of it and build the output list */
    points_number = PyList_Size(points_list); /*  Same as in Python len(_list)  */
    point_dim = PyList_Size(PyList_GetItem(points_list, 0));

    double **points = python_matrix_to_c(points_list, points_number, point_dim);
    double **matrix = jacobi(points, points_number);

    // converting C 2d array to Python list of lists.
    PyObject *py_matrix = PyList_New(0);
    if (py_matrix == NULL)
        other_error();

    for (i = 0; i < points_number + 1; i++)
    {
        PyObject *row_i = PyList_New(0);
        if (row_i == NULL)
            other_error();
        for (j = 0; j < points_number; j++)
        {
            PyObject *float_point = PyFloat_FromDouble((double)matrix[i][j]);
            PyList_Append(row_i, float_point);
        }
        if (PyList_Append(py_matrix, row_i) == -1)
            other_error();
    }

    // free memory

    free_2d(matrix);
    free_2d(points);
    return py_matrix;
}

/*
 * A macro to help us with defining the methods
 */
#define FUNC(_flag, _name, _docstring)                           \
    {                                                            \
#_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) \
    }

static PyMethodDef _methods[] = {
    FUNC(METH_VARARGS, kmeans_fit, "K-means-algorithm"),
    FUNC(METH_VARARGS, wam_fit, "Calculate weighted adjacency matrix"),
    FUNC(METH_VARARGS, ddg_fit, "Calculate diagonal degree matrix"),
    FUNC(METH_VARARGS, lnorm_fit, "Calculate normalized graph laplacian"),
    FUNC(METH_VARARGS, jacobi_fit, "Calculate eiganvalues and eiganvectores"),
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans module",
    NULL,
    -1,
    _methods};

PyMODINIT_FUNC PyInit_spkmeans_module(void)
{
    return PyModule_Create(&_moduledef);
}

// static PyMethodDef _methods[] = {
//     FUNC(METH_VARARGS, fit, "K-means-algorithm"),
//     {NULL, NULL, 0, NULL} /* sint ientinel */
// };

// static struct PyModuleDef _moduledef = {
//     PyModuleDef_HEAD_INIT,
//     "mykmeanssp",
//     NULL,
//     -1,
//     _methods};

// PyMODINIT_FUNC
// PyInit_mykmeanssp(void)
// {
//     return PyModule_Create(&_moduledef);
// }