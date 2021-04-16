from fire import Fire
import numpy as np
import pandas as pd


def build_inputs(fpath):
    with open(fpath, "r") as f:
        n = int(f.readline())
        rows = []
        for line in f.readlines():
            rows.append([int(x) for x in line.split()])
    dist_mat = pd.DataFrame(data=np.array(rows), index=range(n), columns=range(n))
    return n, dist_mat


def upgma(dist_mat, n):
    clusts = set(range(n))
    grph = {i: None for i in range(n)}
    age = np.zeros(n)
    while len(clusts) > 1:
        c_i, c_j = find_closest(dist_mat)


def main(fpath):
    n, dist_mat = build_inputs(fpath)
    tree = upgma()


if __name__ == "__main__":
    Fire(main)
