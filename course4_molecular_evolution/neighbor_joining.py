import numpy as np


def d_star(dmat):
    out = np.zeros(dmat.shape)
    n = dmat.shape[0]
    total_dists = dmat.sum(axis=1)
    for i in range(dmat.shape[0]):
        for j in range(dmat.shape[1]):
            if i == j:
                continue
            out[i, j] = (n-2)*dmat[i, j] - total_dists[i] - total_dists[j]
    return out


def nearest_nebs(d_star):
    return np.argwhere(d_star==np.min(d_star))[0]


def new_dist(dmat, i, j, k):
    assert k != i, 'k = i'
    assert k != j, 'k = j'
    return 0.5 * (dmat[k,i] + dmat[k,j] - dmat[i,j])


def delta_ij(dmat, i, j):
    total_dists = dmat.sum(axis=1)
    return (total_dists[i] - total_dists[j]) / (dmat.shape[0] - 2)


def limb_lengths(dmat, i, j):
    delta = delta_ij(dmat, i, j)
    return 0.5*(dmat[i,j] + delta), 0.5*(dmat[i,j] - delta)
