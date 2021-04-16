from random import choice as rc

from fire import Fire

SAMPLE_ADJ_GRPH = {
    0: [3],
    1: [0],
    2: [1, 6],
    3: [2],
    4: [2],
    5: [4],
    6: [5, 8],
    7: [9],
    8: [7],
    9: [6],
}


def build_adj_grph(grph_file):
    with open(grph_file, "r") as f:
        lines = f.readlines()
    adj_grph = dict()
    for line in lines:
        if "->" not in line:
            continue
        start, ends = line.split("->")
        start = int(start)
        ends = [int(node.strip()) for node in ends.split(",")]
        adj_grph[start] = ends
    return adj_grph


def count_all_edges(adj_grph):
    all_edges = set()
    for start_node, end_nodes in adj_grph.items():
        for end_node in end_nodes:
            all_edges.add((start_node, end_node))
    return all_edges


def find_all_paths(cur_node, adj_grph):
    return set([(cur_node, next_node) for next_node in adj_grph[cur_node]])


def random_walk(start_node, adj_grph, already_walked=None, nodes_unused=None):
    if already_walked is None:
        walk = list()
    else:
        walk = already_walked
    if nodes_unused is None:
        nodes_unused = set()

    cur_node = start_node
    while True:
        options = find_all_paths(cur_node, adj_grph)
        if options.issubset(set(walk)):
            # Nowhere left to go
            return walk, nodes_unused
        valid_options = options - set(walk)
        next_edge = rc(list(valid_options))
        walk.append(next_edge)
        cur_node = next_edge[1]
        has_unused = has_unused_edges(cur_node, walk, adj_grph)
        if has_unused:
            nodes_unused.add(cur_node)
        if not has_unused and cur_node in nodes_unused:
            nodes_unused.remove(cur_node)


def has_unused_edges(node, walk, adj_grph):
    all_paths = find_all_paths(node, adj_grph)
    return len(all_paths - set(walk)) > 1


def walk_previous_cycle(new_start, prev_cycle):
    start_ix = [tup[0] for tup in prev_cycle].index(new_start)
    return prev_cycle[start_ix:] + prev_cycle[:start_ix]


def eulerian_path(adj_grph):
    all_edges = count_all_edges(adj_grph)
    start_node = rc(list(adj_grph.keys()))
    cycle, nodes_unused = random_walk(
        start_node, adj_grph, already_walked=None, nodes_unused=None
    )
    while set(cycle) != all_edges:
        new_start = rc(list(nodes_unused))
        nodes_unused.remove(new_start)
        cycle = walk_previous_cycle(new_start, cycle)
        cycle, nodes_unused = random_walk(
            start_node=new_start,
            adj_grph=adj_grph,
            already_walked=cycle,
            nodes_unused=nodes_unused,
        )
    return cycle


def print_walk(walk):
    assert walk[0][0] == walk[-1][1], "We don't start and end at the same place"
    to_print = [edge[0] for edge in walk]
    to_print.append(walk[0][0])
    to_print = [str(x) for x in to_print]
    print("->".join(to_print))


def main(fpath):
    adj_grph = build_adj_grph(grph_file=fpath)
    cycle = eulerian_path(adj_grph=adj_grph)
    print_walk(cycle)


if __name__ == "__main__":
    Fire(main)
