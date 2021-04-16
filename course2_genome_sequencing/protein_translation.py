from collections import Counter

from fire import Fire
import numpy as np


STOP_CODON = "x"

NUC_2_AMINO = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAU": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGU": "S",
    "AUA": "I",
    "AUC": "I",
    "AUG": "M",
    "AUU": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAU": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAU": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",
    "UAA": STOP_CODON,
    "UAC": "Y",
    "UAG": STOP_CODON,
    "UAU": "Y",
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",
    "UGA": STOP_CODON,
    "UGC": "C",
    "UGG": "W",
    "UGU": "C",
    "UUA": "L",
    "UUC": "F",
    "UUG": "L",
    "UUU": "F",
}

DNA_MAP = {"A": "T", "C": "G", "G": "C", "T": "A"}

AMINO_ACID_MASS = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'I': 113,
    'L': 113,
    'N': 114,
    'D': 115,
    'K': 128,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186,
}


def reverse_complement(dna):
    return ''.join([DNA_MAP[nuc] for nuc in reversed(dna)])


def reading_frames(rna):
    for start in range(3):
        yield rna[start:]


def find_codons(rna):
    codons = [rna[i : i + 3] for i in range(0, len(rna), 3)]
    if len(codons[-1]) != 3:
        codons = codons[:-1]  # strip the incomplete last one
    return codons


def transcribe(rna):
    return rna.replace("T", "U")


def translate(rna, codon_map=NUC_2_AMINO):
    codons = find_codons(rna)
    aminos = [codon_map[codon] for codon in codons]
    if aminos[-1] == "x":
        aminos = aminos[:-1]  # drop stop codon
    return ''.join(aminos)


def find_all_start_ixs(aminos, peptides):
    result = []
    pep_len = len(peptides)
    for poss_start in range(len(aminos) - pep_len + 1):
        if aminos[poss_start : poss_start + pep_len] == peptides:
            result.append(poss_start)
    return result


def find_rna_for_peptides(rna, peptides):
    result = []
    for frame in reading_frames(rna):
        frame = transcribe(frame)
        aminos = translate(frame)
        start_ixs = find_all_start_ixs(aminos, peptides)
        for start_ix in start_ixs:
            rna_start = start_ix * 3
            rna_end = rna_start + (len(peptides) * 3)
            rna_frag = rna[rna_start:rna_end]
            result.append(rna_frag)
    for frame in reading_frames(reverse_complement(rna)):
        frame = transcribe(frame)
        aminos = translate(frame)
        start_ixs = find_all_start_ixs(aminos, peptides)
        for start_ix in start_ixs:
            rna_start = start_ix * 3
            rna_end = rna_start + (len(peptides) * 3)
            result.append(rna[rna_start:rna_end])
    return result


def prefix_mass(peptide, mass_dict=AMINO_ACID_MASS):
    p_mass = np.zeros(len(peptide) + 1, dtype=np.int16)
    for prefix_ix in range(1, len(peptide) + 1):
        sym = peptide[prefix_ix - 1]
        p_mass[prefix_ix] = p_mass[prefix_ix - 1] + mass_dict[sym]
    return p_mass


def linear_spectrum(peptide, mass_dict=AMINO_ACID_MASS):
    p_mass = prefix_mass(peptide, mass_dict)
    linear_spec = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            linear_spec.append(p_mass[j] - p_mass[i])
    return sorted(linear_spec)


def cyclic_spectrum(peptide, mass_dict=AMINO_ACID_MASS):
    p_mass = prefix_mass(peptide, mass_dict)
    peptide_mass = p_mass[-1]
    assert peptide_mass == max(p_mass), "Weirdness here"
    cyclic_spec = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            cut = p_mass[j] - p_mass[i]
            cyclic_spec.append(cut)
            if i > 0 and j < len(peptide):
                cyclic_spec.append(peptide_mass - cut)
    return sorted(cyclic_spec)


def expand(candidates, alphabet):
    expanded = set()
    for cand in candidates:
        for symbol in alphabet:
            expanded.add(cand + symbol)
    return expanded


def mass(peptide, mass_dict=AMINO_ACID_MASS):
    return sum([mass_dict[amino] for amino in peptide])


def parent_mass(spectrum):
    return max(spectrum)


def consistent(peptide, spectrum):
    pep_spec = linear_spectrum(peptide)
    pep_cntr = Counter(pep_spec)
    spec_cntr = Counter(spectrum)
    for wt in pep_cntr.keys():
        if wt not in spec_cntr:
            return False
        elif pep_cntr[wt] > spec_cntr[wt]:
            return False
    return True


def score(peptide, spectrum, seq_type="cyclic"):
    if seq_type == "cyclic":
        pep_spec = cyclic_spectrum(peptide)
    elif seq_type == "linear":
        pep_spec = linear_spectrum(peptide)
    else:
        raise ValueError(f"seq_type can be 'cyclic' OR 'linear', not {seq_type}")
    pep_cntr = Counter(pep_spec)
    spec_cntr = Counter(spectrum)
    score = 0
    for wt in pep_cntr:
        if wt in spec_cntr:
            score += min(pep_cntr[wt], spec_cntr[wt])
    return score


def cyclopep_sequencing(spectrum, mass_dict=AMINO_ACID_MASS):
    cands = set([""])
    final_peps = set()
    while len(cands) > 0:
        cands = expand(candidates=cands, alphabet=mass_dict.keys())
        for peptide in cands.copy():
            if mass(peptide) == parent_mass(spectrum):
                if peptide not in final_peps and cyclic_spectrum(peptide) == spectrum:
                    final_peps.add(peptide)
                cands.remove(peptide)
            elif not consistent(peptide, spectrum):
                cands.remove(peptide)
    return final_peps


def leaderboard_cyclopep_sequencing(spectrum, N, mass_dict=AMINO_ACID_MASS):
    leader_brd = set([""])
    leader = ""
    while len(leader_brd) > 0:
        leader_brd = expand(candidates=leader_brd, alphabet=mass_dict.keys())
        for peptide in leader_brd.copy():
            if mass(peptide) == parent_mass(spectrum):
                if score(peptide, spectrum) > score(leader, spectrum):
                    leader = peptide
            elif mass(peptide) > parent_mass(spectrum):
                leader_brd.remove(peptide)
        leader_brd = trim(leader_brd, spectrum, N)
    return leader


def trim(leader_brd, spectrum, N):
    if len(leader_brd) == 0:
        return leader_brd
    scores = []
    for peptide in leader_brd:
        pep_score = score(peptide, spectrum, seq_type="linear")
        scores.append((peptide, pep_score))
    scores.sort(key=lambda tup: tup[1], reverse=True)
    N = min(N, len(leader_brd)-1)
    cutoff = scores[N][1]
    return set([tup[0] for tup in scores if tup[1] >= cutoff])


def format_peptides_as_weights(peps, mass_dict=AMINO_ACID_MASS):
    result = set()
    for pep in peps:
        result.add("-".join([str(mass_dict[amino]) for amino in pep]))
    return result


def spectral_convolution(spectrum):
    result = []
    spectrum.sort(reverse=True)
    start_ix = 0
    while start_ix < len(spectrum) - 1:
        for next_ix in range(start_ix+1, len(spectrum)):
            result.append(spectrum[start_ix] - spectrum[next_ix])
        start_ix += 1
    return [str(num) for num in result if num!= 0]


def main(spectrum):
    spectrum = [int(x) for x in spectrum.split()]
    print(' '.join(spectral_convolution(spectrum)))


if __name__ == "__main__":
    Fire(main)
