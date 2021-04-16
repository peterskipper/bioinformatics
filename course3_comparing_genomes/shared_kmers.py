from collections import defaultdict
from itertools import islice

DNA_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def window(seq, n):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def reverse_complement(dna):
    return ''.join([DNA_MAP[nuc] for nuc in reversed(dna)])


def shared_kmers(strand1, strand2, k):
    seen_kmers = defaultdict(list)
    for ix, kmer in enumerate(window(strand1, k)):
        seen_kmers[''.join(kmer)].append(ix)
        seen_kmers[reverse_complement(kmer)].append(ix)

    shared_kmers = 0
    for kmer in window(strand2, k):
        if ''.join(kmer) in seen_kmers:
            shared_kmers += 1

    return shared_kmers, seen_kmers
