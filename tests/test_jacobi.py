import sys

from tests.load_module import load_test_data
sys.path.append("../Spectral-Clustering-Algorithm/src/")
import spkmeans_module as spkm
import pytest


def format_arr_to_4digits(arr):
    for row in arr:
        for i in range(len(row)):
            row[i] = float(f'{row[i]:.4f}')
    return arr


test_data = load_test_data("jacobi")


@pytest.mark.parametrize("input,correct_output", test_data)
def test_jacobi(input, correct_output):
    output = format_arr_to_4digits(spkm.jacobi_fit(input))
    # check without order!!! use http://www.math.utoledo.edu/~codenth/Linear_Algebra/Calculators/jacobi_algorithm.html
    assert (output == correct_output)
