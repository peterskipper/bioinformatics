TEST_STR = "panamabananas$"


def make_suffixes(string):
    result = []
    for ix in range(len(string)):
        result.append((string[ix:], ix))
    return result


def suffix_arr(string):
    suffix_tuples = make_suffixes(string)
    suffix_tuples.sort(key=lambda tup: tup[0])
    return [tup[1] for tup in suffix_tuples]


def cyclic_rots(string):
    result = []
    for ix in range(len(string)):
        result.append(string[ix:] + string[:ix])
    return result


def bwt(string):
    cycles = cyclic_rots(string)
    cycles.sort()
    return [cycle[-1] for cycle in cycles]


def enumerate_col(string):
    seen = {}
    result = []
    for char in string:
        if char in seen:
            seen[char].append(char)
        else:
            seen[char] = [char]
        result.append(f'{char}{len(seen[char])}')
    return result


def reconstruct_path(first_col, last_col):
    assert first_col[0] == "$1", "oops"
    result = []
    first_tok = first_col[0]
    while len(result) != len(first_col):
        cur_ix = last_col.index(first_tok)
        result.append(first_col[cur_ix])
        first_tok = first_col[cur_ix]
    return result


def invert_bwt(bwt_str):
    first_col = ''.join(sorted([tok for tok in bwt_str]))
    first_col_tokens = enumerate_col(first_col)
    last_col_tokens = enumerate_col(bwt_str)
    path = reconstruct_path(first_col_tokens, last_col_tokens)
    return ''.join([tok[0] for tok in path])

