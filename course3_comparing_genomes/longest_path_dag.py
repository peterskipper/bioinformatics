from fire import Fire
import numpy as np
from queue import SimpleQueue


def build_graph(lines):
    grph = dict()
    for line in lines:
        start, end = line.split("->")
        end, wt = end.split(":")
        wt = int(wt)
        if start in grph:
            grph[start].append((end, wt))
        else:
            grph[start] = [(end, wt)]
    return grph


def find_all_nodes(grph):
    all_nodes = set(grph.keys())
    for node in grph:
        for neigh, _ in grph[node]:
            if neigh not in all_nodes:
                all_nodes.add(neigh)
    return list(all_nodes)


def add_nodes_to_graph(nodes, grph):
    """Adds all nodes to grph, whether or not they have outgoing edges"""
    for node in nodes:
        if node not in grph:
            grph[node] = []
    return grph


def is_connected(grph, start, end):
    if start == end:
        return True

    visited = set()
    visited.add(start)

    q = SimpleQueue()
    q.put(start)

    while not q.empty():
        node = q.get()
        for neigh, _ in grph[node]:
            if neigh == end:
                return True
            if neigh not in visited:
                q.put(neigh)
                visited.add(neigh)
    return False


def prune(grph, nodes, source=None, sink=None):
    if source is None and sink is None:
        raise AssertionError("Only 1 of source/sink can be None")
    elif source is not None and sink is not None:
        raise AssertionError("Only 1 of source/sink can be None")

    if source is not None:
        return [node for node in nodes if is_connected(grph, source, node)]
    if sink is not None:
        return [node for node in nodes if is_connected(grph, node, sink)]


def find_in_degrees(grph, nodes):
    in_degrees = {node: 0 for node in nodes}
    for start in grph.keys():
        for end, wt in grph[start]:
            if end in in_degrees:
                in_degrees[end] += 1
    return in_degrees


def kahn_topo_sort(grph, nodes, source=None, sink=None):
    result = []
    in_degrees = find_in_degrees(grph, nodes)
    no_inpaths = set([node for node, indg in in_degrees.items() if indg == 0])
    q = SimpleQueue()
    if source is not None:
        # Make sure source is at the start of topo sort
        q.put(source)
        no_inpaths = no_inpaths - set([source])
    for node in no_inpaths:
        q.put(node)
    while not q.empty():
        node = q.get()
        result.append(node)
        for neigh, wt in grph[node]:
            if neigh not in in_degrees:
                continue
            in_degrees[neigh] -= 1
            if in_degrees[neigh] == 0:
                q.put(neigh)
    if sink is not None and result[-1] != sink:
        # if sink is not at the end of topo sort, place it there
        sink_index = result.index(sink)
        result.append(sink)
        result = result[:sink_index] + result[sink_index + 1 :]
    return result


def find_predecessors(grph, topo_sort, ix):
    cur_node = topo_sort[ix]
    result = []
    for node in topo_sort[:ix]:
        direct_connect = set([n for n, wt in grph[node]])
        if cur_node in direct_connect:
            result.append(node)
    return result


def find_weight(grph, start, end):
    connected = grph[start]
    for node, wt in connected:
        if node == end:
            return wt
    raise AssertionError("something is amiss")


def find_longest_path(source, sink, grph):
    all_nodes = find_all_nodes(grph)
    grph = add_nodes_to_graph(all_nodes, grph)
    valid_nodes = prune(grph=grph, nodes=all_nodes, source=source)
    valid_nodes = prune(grph=grph, nodes=valid_nodes, sink=sink)
    topo_sorted_nodes = kahn_topo_sort(
        grph=grph, nodes=valid_nodes, source=source, sink=sink
    )

    assert topo_sorted_nodes[0] == source, "shenanigans"
    path_length = {}
    longest_path = {}
    for ix, node in enumerate(topo_sorted_nodes):
        if ix == 0:
            path_length[node] = 0
            longest_path[node] = [node]
            continue
        preds = find_predecessors(grph, topo_sorted_nodes, ix)
        longest = -np.inf
        for pred in preds:
            pth_length = path_length[pred] + find_weight(grph, pred, node)
            if pth_length > longest:
                longest = pth_length
                long_pth = longest_path[pred][:] + [node]
        path_length[node] = longest
        longest_path[node] = long_pth
    return path_length[sink], longest_path[sink]


def main(fpath):
    with open(fpath, 'r') as f:
        lines = f.readlines()
    source = lines[0].strip()
    sink = lines[1].strip()
    grph = build_graph(lines[2:])

    all_nodes = find_all_nodes(grph)
    pth_len, pth = find_longest_path(source, sink, grph)
    print(pth_len)
    pth = [str(x) for x in pth]
    print("->".join(pth))


if __name__ == "__main__":
    Fire(main)
