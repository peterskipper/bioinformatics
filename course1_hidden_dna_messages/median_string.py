from itertools import product

from fire import Fire
import numpy as np

from utils import window
from hamming import hamming_dist


def min_dist_strand(pattern, strand):
    min_dist = np.inf
    n = len(pattern)
    for wind in window(strand, n):
        dist = hamming_dist(pattern, wind)
        if dist < min_dist:
            min_dist = dist
    return min_dist


def min_dist(pattern, dna):
    result = 0
    for strand in dna:
        result += min_dist_strand(pattern, strand)
    return result


def median_string_brute(dna_file, k):
    with open(dna_file, 'r') as f:
        dna = f.read().splitlines()
    best_distance = np.inf
    median = []
    for median_cand in product(['A','C','G','T'], repeat=k):
        median_cand = ''.join(median_cand)
        if len(set(median_cand)) == 1:
            print(f"On {median_cand}")
        d = min_dist(median_cand, dna)
        if d < best_distance:
            best_distance = d
            median = [median_cand]
        elif d == best_distance:
            median.append(median_cand)
    return median


if __name__ == "__main__":
    Fire(median_string_brute)

