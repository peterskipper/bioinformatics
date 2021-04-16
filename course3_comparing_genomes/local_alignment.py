from fire import Fire
import numpy as np
import pandas as pd

import global_alignment as ga


def map_backtrack_ix_local(bt_ix):
    if bt_ix == 0:
        return "source"
    elif bt_ix == 1:
        return "down"
    elif bt_ix == 2:
        return "right"
    elif bt_ix == 3:
        return "diag"


def local_lcs(first, second, score_mat, indel_pen):
    score = np.zeros((len(first) + 1, len(second) + 1))
    backtrack = np.empty((len(first) + 1, len(second) + 1), dtype="object")
    backtrack[0][0] = ""
    for i in range(1, len(first) + 1):
        backtrack[i, 0] = "source"
    for j in range(1, len(second) + 1):
        backtrack[0, j] = "source"
    for i in range(1, len(first) + 1):
        for j in range(1, len(second) + 1):
            amino1 = first[i - 1]
            amino2 = second[j - 1]
            scores = np.array(
                [
                    0,  # free taxi ride
                    score[i - 1, j] + indel_pen,
                    score[i, j - 1] + indel_pen,
                    score[i - 1, j - 1] + score_mat.loc[amino1, amino2],
                ]
            )
            score[i, j] = max(scores)
            bt_ix = np.argmax(scores)
            backtrack[i, j] = map_backtrack_ix_local(bt_ix)
    return score, backtrack


def do_backtrack(row_ix, col_ix, backtrack, first, second):
    first_align = []
    second_align = []
    cur = backtrack[row_ix, col_ix]
    while cur != "source":
        if cur == "down":
            first_align.append(first[row_ix - 1])
            second_align.append("-")
            row_ix -= 1
        elif cur == "right":
            first_align.append("-")
            second_align.append(second[col_ix - 1])
            col_ix -= 1
        elif cur == "diag":
            first_align.append(first[row_ix - 1])
            second_align.append(second[col_ix - 1])
            row_ix -= 1
            col_ix -= 1
        else:
            raise AssertionError(f"cur is {cur}")
        cur = backtrack[row_ix, col_ix]
    first_align.reverse()
    second_align.reverse()
    return first_align, second_align


def score_alignments(first_align, second_align, match=1, mismatch=-1, indel_pen=-2):
    assert len(first_align) == len(second_align), "alignments wrong length"
    score = 0
    for first, second in zip(first_align, second_align):
        if first == "-" or second == "-":
            score += indel_pen
        elif first == second:
            score += match
        else:
            score += mismatch
    return score


def main(first, second, indel_pen=-5):
    pam250 = ga.load_score_matrix("pam_250_score_matrix.txt")
    score, backtrack = local_lcs(first, second, pam250, indel_pen)
    highest_score = int(np.max(score))
    row_start_ix = np.argmax(np.max(score, axis=1))
    col_start_ix = np.argmax(np.max(score, axis=0))
    first_align, second_align = do_backtrack(
        row_start_ix, col_start_ix, backtrack, first, second
    )
    print(highest_score)
    print("".join(first_align))
    print("".join(second_align))


if __name__ == "__main__":
    Fire(main)
