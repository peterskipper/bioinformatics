from fire import Fire
import numpy as np
import pandas as pd


def load_score_matrix(score_matrix_pth):
    with open(score_matrix_pth, "r") as f:
        cols = f.readline().split()
        index = []
        data = []
        for line in f.readlines():
            tokens = line.split()
            index.append(tokens[0])
            data.append([int(x) for x in tokens[1:]])
        return pd.DataFrame(data=data, index=index, columns=cols)


def map_backtrack_ix(bt_ix):
    if bt_ix == 0:
        return "down"
    elif bt_ix == 1:
        return "right"
    elif bt_ix == 2:
        return "diag"
    else:
        raise AssertionError("just...why")


def lcs(first, second, score_mat, indel_pen):
    score = np.zeros((len(first) + 1, len(second) + 1))
    backtrack = np.empty((len(first) + 1, len(second) + 1), dtype="object")
    backtrack[:] = ""
    for i in range(1, len(first) + 1):
        score[i, 0] = score[i - 1, 0] + indel_pen
        backtrack[i, 0] = "down"
    for j in range(1, len(second) + 1):
        score[0, j] = score[0, j - 1] + indel_pen
        backtrack[0, j] = "right"
    for i in range(1, len(first) + 1):
        for j in range(1, len(second) + 1):
            amino1 = first[i - 1]
            amino2 = second[j - 1]
            scores = np.array(
                [
                    score[i - 1, j] + indel_pen,
                    score[i, j - 1] + indel_pen,
                    score[i - 1, j - 1] + score_mat.loc[amino1, amino2],
                ]
            )
            score[i, j] = max(scores)
            bt_ix = np.argmax(scores)
            backtrack[i, j] = map_backtrack_ix(bt_ix)
    return score, backtrack


def build_lcs_directions(backtrack):
    row_ix, col_ix = [val - 1 for val in backtrack.shape]
    result = []
    cur = backtrack[row_ix, col_ix]
    while cur != "":
        result.append(cur)
        if cur == "diag":
            row_ix -= 1
            col_ix -= 1
        elif cur == "down":
            row_ix -= 1
        elif cur == "right":
            col_ix -= 1
        else:
            raise AssertionError("How can this be")
        cur = backtrack[row_ix, col_ix]
    result.reverse()
    return result


def write_lcs(first, second, lcs_directions):
    first_ix = 0
    second_ix = 0
    first_align = []
    second_align = []
    for direc in lcs_directions:
        if direc == "down":
            first_align.append(first[first_ix])
            first_ix += 1
            second_align.append("-")
        elif direc == "right":
            first_align.append("-")
            second_align.append(second[second_ix])
            second_ix += 1
        elif direc == "diag":
            first_align.append(first[first_ix])
            first_ix += 1
            second_align.append(second[second_ix])
            second_ix += 1
    assert first_ix == len(first), "First is messed up"
    assert second_ix == len(second), "Second is messed up"
    return first_align, second_align


def main(first, second, indel_pen=-5):
    blosum = load_blosum("blosum_score_matrix.txt")
    score, backtrack = lcs(first, second, blosum, indel_pen)
    max_score = int(score[-1][-1])
    lcs_directions = build_lcs_directions(backtrack)
    first_align, second_align = write_lcs(first, second, lcs_directions)
    print("Max Score: ", max_score)
    print("First Align: ", "".join(first_align))
    print("Second Align: ", "".join(second_align))


if __name__ == "__main__":
    Fire(main)
