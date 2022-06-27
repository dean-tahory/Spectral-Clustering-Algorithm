#include "spkmeans.c"
#include <Python.h>

// method to convert python list of lists to 2d array in C
static double **python_2d_array_to_c(PyObject *_list)
{
    PyObject *point, *point_value;
    Py_ssize_t i, j, n, point_size;

    /* check if this list */
    if (!PyList_Check(_list))
        return NULL;
    /* Get the size of it and build the output list */
    n = PyList_Size(_list); /*  Same as in Python len(_list)  */

    // creating the points array where each point is array of double.
    double **points = calloc(n, sizeof(double *));
    for (i = 0; i < n; i++)
    {
        point = PyList_GetItem(_list, i);
        if (!PyList_Check(point))
        { /* We only adding lists of points */
            continue;
        }
        point_size = PyList_Size(point);
        points[i] = malloc(sizeof points[i] * point_size);
        /* Check that we got the memory from OS. In the assert - a string has a true value */
        assert(points[i] != NULL && "Problem in python_2d_array_to_c()");
        for (j = 0; j < point_size; j++)
        {
            point_value = PyList_GetItem(point, j);
            if (!PyFloat_Check(point_value))
                continue; /* Skip non-float values */

            points[i][j] = PyFloat_AS_DOUBLE(point_value); /* Convert a Python float object to double */
            if (points[i][j] == -1 && PyErr_Occurred())
            {
                /* Flaot too big to fit in a C double, bail out */
                free(points[i]);
                other_error();
            }
        }
    }
    return points;
}

// our main module API method - getting args from python object, converting them to C objects and work on them.
static PyObject *fit(PyObject *self, PyObject *args)
{
    int K, max_iter, points_length, point_length, i;
    double eps;
    PyObject *points_list, *centroids_list;

    // binding the variables we declared to the args we get from the python module.
    if (!PyArg_ParseTuple(args, "iiOOiid:fit", &K, &max_iter, &points_list, &centroids_list, &points_length, &point_length, &eps))
    {
        return NULL;
    }

    double **points = python_2d_array_to_c(points_list);
    double **initial_centroids = python_2d_array_to_c(centroids_list);

    // calling the kmeans C algorithm (new version that also gets initial centroids from the python module)
    double **centroids = k_means(K, max_iter, points, initial_centroids, points_length, point_length, eps);

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
        for (j = 0; j < point_length; j++)
        {
            PyObject *float_point = PyFloat_FromDouble(centroids[i][j]);
            PyList_Append(new_centroid_point, float_point);
        }
        if (PyList_Append(new_centroids, new_centroid_point) == -1)
            other_error();
    }
    return new_centroids;
}

/*
 * A macro to help us with defining the methods
 */
#define FUNC(_flag, _name, _docstring)                           \
    {                                                            \
#_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) \
    }

static PyMethodDef _methods[] = {
    FUNC(METH_VARARGS, spk, "Spectral-Clustering-Algorithm"),
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _methods};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    return PyModule_Create(&_moduledef);
}