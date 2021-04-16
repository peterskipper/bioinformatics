from fire import Fire
from numpy import argmax


def max_skew_ix(nucleotides):
    cur_skew = 0
    skew_arr = []
    for ix, nuc in enumerate(nucleotides):
        skew_arr.append(cur_skew)
        if nuc.lower() == "g":
            cur_skew += 1
        if nuc.lower() == "c":
            cur_skew -= 1
    skew_arr.append(cur_skew)
    return argmax(skew_arr)


if __name__ == "__main__":
    Fire(max_skew_ix)
