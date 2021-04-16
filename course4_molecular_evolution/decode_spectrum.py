import networkx as nx
import numpy as np


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
    # Fake Peptides for toy code examples
    'X': 4,
    'Z': 5
}

TEST_INPUT = [57, 71, 154, 185, 301, 332, 415, 429, 486]


def build_nodes(input_lst):
    return [0] + sorted(input_lst)


def build_mass2amino():
    return {v: k for k, v in AMINO_ACID_MASS.items()}


def build_graph(nodes, mass2amino):
    G = nx.DiGraph()
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            diff = nodes[j] -nodes[i]
            if diff in mass2amino:
                G.add_edge(nodes[i], nodes[j], label=mass2amino[diff])
    return G


def format_graph(G):
    edge_lst = sorted(G.edges, key=lambda tup: tup[0])
    for source, dest in edge_lst:
        lbl = G.get_edge_data(source, dest)['label']
        print(f"{source}->{dest}:{lbl}")


def print_graph_from_spectrum(spectrum):
    nodes = build_nodes(spectrum)
    mass2amino = build_mass2amino()
    G = build_graph(nodes, mass2amino)
    format_graph(G)


def weigh_peptide(peptide):
    if peptide == "":
        return 0
    wt = 0
    for pep in peptide:
        wt += AMINO_ACID_MASS[pep]
    return wt


def ideal_spectrum(peptide):
    spectrum = []
    for cut in range(len(peptide)):
        prefix = peptide[:cut]
        suffix = peptide[cut:]
        spectrum.append(weigh_peptide(prefix))
        spectrum.append(weigh_peptide(suffix))
    return sorted(spectrum)


def path_to_peptide(G, path):
    return ''.join([G.get_edge_data(source, dest)['label'] for source, dest in path])


def decode_ideal_spectrum(spectrum):
    nodes = build_nodes(spectrum)
    mass2amino = build_mass2amino()
    G = build_graph(nodes, mass2amino)
    paths = nx.all_simple_paths(G, source=nodes[0], target=nodes[-1])
    for path in map(nx.utils.pairwise, paths):
        peptide = path_to_peptide(G, path)
        if ideal_spectrum(peptide) == nodes:
            return peptide


def peptide2vector(peptide):
    prefix_wts = []
    for cut in range(1, len(peptide)+1):
        prefix = peptide[:cut]
        prefix_wts.append(weigh_peptide(prefix))
    prefix_wts = sorted(prefix_wts)
    pep_vector = np.zeros(prefix_wts[-1])
    for pw in prefix_wts:
        pep_vector[pw-1] = 1
    return pep_vector
