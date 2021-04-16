from fire import Fire
from itertools import combinations
import numpy as np


def limb_length(j, dist_mat):
    valid_ixs = list(set(range(dist_mat.shape[0])) - set([j]))
    results = []
    for (i, k) in combinations(valid_ixs, 2):
        res = (dist_mat[i, j] + dist_mat[j, k] - dist_mat[i, k]) / 2
        results.append(res)
    return int(min(results))


def build_inputs(fpath):
    with open(fpath, "r") as f:
        n = int(f.readline())
        j = int(f.readline())
        rows = []
        for line in f.readlines():
            rows.append([int(x) for x in line.split()])
    dist_mat = np.array(rows)
    return n, j, dist_mat


def main(fpath):
    n, j, dist_mat = build_inputs(fpath)
    return limb_length(j, dist_mat)


if __name__ == "__main__":
    Fire(main)
