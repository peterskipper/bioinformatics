from utils import window


def how_likely_is_motif(profile, motif):
    assert profile.shape[1] == len(motif), "motif length and profile cols don't match"
    prob = 1
    for ix in range(len(motif)):
        nuc = motif[ix]
        prob *= profile.loc[nuc, ix]
    return prob


def find_most_probable_motif(profile, strand, k):
    highest_prob = -1
    result = None
    for motif in window(strand, n=k):
        prob = how_likely_is_motif(profile, motif)
        if prob > highest_prob:
            highest_prob = prob
            result = motif
    return motif


def motifs(profile, dna, k):
    motifs = []
    for strand in dna:
        most_prob = find_most_probable_motif(profile, strand, k)
        motifs.append(most_prob)
    return motifs


def score(motifs):
    consensus = find_consensus(motifs)
    score = 0
    for ix in range(len(motifs[0])):
        target = 
