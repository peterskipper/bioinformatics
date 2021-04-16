import numpy as np

DATA = [(2,8), (2,5), (6,9), (7,5), (5,2)]
CENTERS = [(3,5), (5,4)]
HIDDEN_MAT = [[0.5, 0.3, 0.8, 0.4, 0.9],
              [0.5, 0.7, 0.2, 0.6, 0.1]]


def d(data_pt, center):
    return np.linalg.norm(data_pt-center)


def newton_inv(data, centers):
    hidden_mat = np.zeros((len(centers), len(data)))
    for i, center in enumerate(centers):
        for j, datum in enumerate(data):
            inv_dist = 1/(d(datum, center)**2)
            all_inv_dists = sum([1/(d(datum, center)**2) for center in centers])
            hidden_mat[i, j] = inv_dist/all_inv_dists
    return hidden_mat


def weighted_center_grav(hidden_mat, data):
    all_results = []
    result = None
    for i in range(hidden_mat.shape[0]):
        if result is not None:
            all_results.append(result)
        result = []
        for j in range(len(data[0])):
            jth = [datum[j] for datum in data]
            weighted = np.dot(hidden_mat[i,:], jth)
            result.append(weighted / sum(hidden_mat[i,:]))
    all_results.append(result)
    return all_results


if __name__ == "__main__":
    DATA = [np.array(dt) for dt in DATA]
    CENTERS = [np.array(center) for center in CENTERS]
    HIDDEN_MAT = np.array(HIDDEN_MAT)
    print(weighted_center_grav(HIDDEN_MAT, DATA))

