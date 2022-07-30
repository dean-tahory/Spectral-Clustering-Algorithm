import json


def load_test_data(test_type: str):
    try:
        f = open("../Spectral-Clustering-Algorithm/tests/test_data.json", 'r')
        raw_data = json.load(f)
        data = [(raw_data[case]["input"], raw_data[case][test_type + "_output"]) for case in raw_data.keys()]
        return data
    except FileNotFoundError as fnf_error:
        print(fnf_error)
