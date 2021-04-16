import numpy as np

DATA = [(1,6), (5,6), (8,7), (1,3), (3,4), (5,2), (7,1), (10,3)]
CENTERS = [(2,4), (6,7), (7,3)]

def d(data_pt, centers):
    dists = [np.linalg.norm(data_pt-center) for center in centers]
    return min(dists)


def max_dist(data, centers):
    data = [np.array(dt) for dt in data]
    centers = [np.array(center) for center in centers]
    all_dists = [d(data_pt, centers) for data_pt in data]
    return max(all_dists)


def transform_array(lst):
    return [np.array(x) for x in lst]


def distortion(data, centers):
    data = transform_array(data)
    centers = transform_array(centers)
    all_dists = [d(data_pt, centers) for data_pt in data]
    squared = np.square(all_dists)
    return np.mean(squared)

