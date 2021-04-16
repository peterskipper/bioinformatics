from fire import Fire

TEST_PERM = [-3, 4, 1, 5, -2]


def perm_from_string(string):
    return [int(x) for x in string.split()]


def make_string(elem):
    if elem > 0:
        return f"+{elem}"
    else:
        return f"{elem}"


def print_perm(perm):
    print(" ".join([make_string(elem) for elem in perm]))


def make_reversal(perm, start_ix, end_ix):
    start = perm[:start_ix]
    end = perm[end_ix:]
    rev = reversed(perm[start_ix:end_ix])
    rev = [elem * -1 for elem in rev]
    return start + rev + end


def greedy_reversal(perm, verbose=False):
    greedy_rev_dist = 0
    for elem in range(1, len(perm) + 1):
        cur_ix = elem - 1
        if perm[cur_ix] != elem:
            try:
                end_ix = perm.index(elem) + 1
            except ValueError:
                end_ix = perm.index(elem * -1) + 1
            perm = make_reversal(perm, start_ix=cur_ix, end_ix=end_ix)
            greedy_rev_dist += 1
            if verbose:
                print_perm(perm)
        if perm[cur_ix] < 0:
            perm[cur_ix] *= -1
            greedy_rev_dist += 1
            if verbose:
                print_perm(perm)
    return greedy_rev_dist


def count_breakpoints(perm):
    breaks = 0
    possible_breaks = [0] + perm + [len(perm) + 1]
    for nbr1, nbr2 in zip(possible_breaks[:-1], possible_breaks[1:]):
        if nbr2 - nbr1 != 1:
            breaks += 1
    return breaks
