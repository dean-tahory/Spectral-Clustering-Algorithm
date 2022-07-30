import sys
sys.path.append("../Spectral-Clustering-Algorithm/src/")
import spkmeans_module as spkm
import pytest

# TODO Read Clean code best tips and tricks in code summaries folder
# TODO learn how to create unit test for this project - https://shorturl.at/aBFX4


def format_arr_to_4digits(arr):
    for row in arr:
        for i in range(len(row)):
            row[i] = float(f'{row[i]:.4f}')
    return arr


def print_2d_list(arr):
    for row in arr:
        for i in range(len(row)):
            row[i] = f'{row[i]:.4f}'
        print(','.join(row))


input = [[1.0, 2.0, 3.0],
         [2.0, 4.0, 6.0],
         [3.0, 6.0, 8.0]]
correct_output = [
    [0.0000, 0.1540, 0.0349],
    [0.1540, 0.0000, 0.2231],
    [0.0349, 0.2231, 0.0000]]

check_list = [(input, correct_output)]


@pytest.mark.parametrize("input,correct_output", check_list)
def test_wam(input, correct_output):
    output = format_arr_to_4digits(spkm.wam_fit(input))
    assert (output == correct_output), f"should be {correct_output}"
    print("test #1 success")


def test_wam2():
    assert 3 < 44

    test_wam(input, correct_output)
