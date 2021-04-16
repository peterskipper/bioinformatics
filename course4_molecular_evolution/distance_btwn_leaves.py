from fire import Fire
import networkx as nx
import numpy as np


def make_edge_w_weights(line):
    start, end_chunk = line.split("->", maxsplit=1)
    end, weight = end_chunk.split(":", maxsplit=1)
    return int(start), int(end), int(weight)


def build_graph(lines):
    weighted_edges = [make_edge_w_weights(line) for line in lines]
    G = nx.Graph()
    G.add_weighted_edges_from(weighted_edges)
    return G


def find_leaves(G):
    return sorted([node for node in G.nodes if G.degree[node] == 1])


def make_distance_matrix(G, num_leaves):
    dist_mat = np.zeros((num_leaves, num_leaves))
    lengths = nx.all_pairs_dijkstra_path_length(G)
    for start, length_dct in lengths:
        if start >= num_leaves:
            continue
        for end in range(num_leaves):
            dist_mat[start, end] = int(length_dct[end])
    return dist_mat


def format_dist_mat(dist_mat):
    for row in range(len(dist_mat)):
        print("\t".join([str(int(d)) for d in dist_mat[row]]))


def main(fpath):
    with open(fpath, "r") as f:
        num_leaves = int(f.readline())
        lines = [line.strip() for line in f.readlines()]
    G = build_graph(lines)
    dist_mat = make_distance_matrix(G, num_leaves)
    format_dist_mat(dist_mat)


if __name__ == "__main__":
    Fire(main)
