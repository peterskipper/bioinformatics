from fire import Fire


def manhattan(n, m, down, right):
    longest = dict()
    longest[(0, 0)] = 0
    for i in range(1, n + 1):
        longest[(i, 0)] = longest[(i - 1, 0)] + down[i - 1][0]
    for j in range(1, m + 1):
        longest[(0, j)] = longest[(0, j - 1)] + right[0][j - 1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            longest[(i, j)] = max(
                longest[(i - 1, j)] + down[i - 1][j],
                longest[(i, j - 1)] + right[i][j - 1],
            )
    return longest[n, m]


def parse_vars(fpath):
    with open(fpath, 'r') as f:
        n, m = [int(val) for val in f.readline().split()]
        down = []
        for _ in range(n):
            line = f.readline()
            row = [int(x) for x in line.split()]
            down.append(row)
        assert f.readline() == "-\n", "off kilter"
        right = []
        for _ in range(n + 1):
            line = f.readline()
            row = [int(x) for x in line.split()]
            right.append(row)
    return n, m, down, right


def main(fpath):
    n, m, down, right = parse_vars(fpath)
    return manhattan(n, m, down, right)


if __name__ == "__main__":
    Fire(main)
