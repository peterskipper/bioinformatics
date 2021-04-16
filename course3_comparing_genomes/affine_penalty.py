import re

from fire import Fire
import numpy as np
import pandas as pd

import global_alignment as ga


def alignment_affine_score(first_align, second_align, match=1, mismatch=-1, gap_open=-4, gap_extend=-1):
    assert len(first_align) == len(second_align), "bloody 'ell"
    score = 0
    gap = False
    score += score_gaps(first_align, gap_open, gap_extend)
    score += score_gaps(second_align, gap_open, gap_extend)
    for f, s in zip(first_align, second_align):
        if f == "-" or s == "-":
            continue
        elif f == s:
            score += match
        else:
            score += mismatch
    return score


def score_gaps(align, gap_open, gap_extend):
    score = 0
    for gap in re.findall("\-+", align):
        score += score_gap(gap, gap_open, gap_extend)
    return score


def score_gap(gap, gap_open, gap_extend):
    return gap_open + (len(gap) - 1) * gap_extend


def init_score(first, second):
    return np.zeros((len(first) + 1, len(second) + 1))


def init_backtrack(first, second):
    backtrack = np.empty((len(first) + 1, len(second) + 1), dtype="object")
    backtrack[:] = ""
    backtrack[0][0] = "source"
    return backtrack


def build_lower_scores(lower, middle, i, j, sigma, epsilon):
    return np.array([lower[i - 1, j] - epsilon, middle[i - 1, j] - sigma])


def build_upper_scores(upper, middle, i, j, sigma, epsilon):
    return np.array([upper[i, j - 1] - epsilon, middle[i, j - 1] - sigma])


def build_middle_scores(lower, middle, upper, i, j, score_mat, amino1, amino2):
    return np.array([lower[i, j], middle[i - 1, j - 1] + score_mat.loc[amino1, amino2], upper[i, j]])


def map_bt_lower(lower_scores):
    bt_ix = np.argmax(lower_scores)
    if bt_ix == 0:
        return "down"
    elif bt_ix == 1:
        return "go_lower_middle"


def map_bt_upper(upper_scores):
    bt_ix = np.argmax(upper_scores)
    if bt_ix == 0:
        return "right"
    elif bt_ix == 1:
        return "go_upper_middle"


def map_bt_middle(middle_scores):
    bt_ix = np.argmax(middle_scores)
    if bt_ix == 0:
        return "go_middle_lower"
    elif bt_ix == 1:
        return "diag"
    elif bt_ix == 2:
        return "go_middle_upper"


def affine_lcs(first, second, score_mat, sigma, epsilon):
    score_lower = init_score(first, second)
    score_lower[0, :] = -np.inf
    score_middle = init_score(first, second)
    score_upper = init_score(first, second)
    score_upper[:, 0] = -np.inf

    bt_lower = init_backtrack(first, second)
    bt_middle = init_backtrack(first, second)
    bt_upper = init_backtrack(first, second)

    for i in range(1, len(first) + 1):
        for j in range(1, len(second) + 1):
            amino1 = first[i - 1]
            amino2 = second[j - 1]

            lower_scores = build_lower_scores(score_lower, score_middle, i, j, sigma, epsilon)
            upper_scores = build_upper_scores(score_upper, score_middle, i, j, sigma, epsilon)
            middle_scores = build_middle_scores(
                score_lower, score_middle, score_upper, i, j, score_mat, amino1, amino2
            )

            score_lower[i, j] = max(lower_scores)
            score_upper[i, j] = max(upper_scores)
            score_middle[i, j] = max(middle_scores)

            bt_lower[i, j] = map_bt_lower(lower_scores)
            bt_upper[i, j] = map_bt_upper(upper_scores)
            bt_middle[i, j] = map_bt_middle(middle_scores)
    return score_lower, score_middle, score_upper, bt_lower, bt_middle, bt_upper


def find_score_and_matrix_to_start_backtrack(score_lower, score_upper, score_middle):
    final_scores = np.array([score_lower[-1][-1], score_upper[-1][-1], score_middle[-1][-1]])
    highest_score = np.max(final_scores)
    highest_ix = np.argmax(final_scores)
    if highest_ix == 0:
        return highest_score, "lower"
    elif highest_ix == 1:
        return highest_score, "upper"
    else:
        return highest_score, "middle"


def do_backtrack(bt_dict, start_key, first, second):
    first_align = []
    second_align = []
    cur_matrix = bt_dict[start_key]
    cur_key = start_key
    row_ix = cur_matrix.shape[0] - 1
    col_ix = cur_matrix.shape[1] - 1
    cur = cur_matrix[row_ix][col_ix]
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
        elif cur == "go_lower_middle":
            cur_key = "lower"
        elif cur == "go_upper_middle":
            cur_key = "upper"
        elif cur == "go_middle_lower":
            cur_key = "middle"
            row_ix -= 1
        elif cur == "go_middle_upper":
            cur_key = "middle"
            col_ix -= 1
        else:
            raise AssertionError(
                f"cur is {cur}\ncur_key is {cur_key}\nrow_ix is {row_ix}\ncol_ix is {col_ix}"
            )
        cur_matrix = bt_dict[cur_key]
        cur = cur_matrix[row_ix, col_ix]
    first_align.reverse()
    second_align.reverse()
    return first_align, second_align


def main(first, second, sigma=11, epsilon=1):
    blosum = ga.load_score_matrix("blosum_score_matrix.txt")
    score_lower, score_middle, score_upper, bt_lower, bt_middle, bt_upper = affine_lcs(
        first, second, blosum, sigma, epsilon
    )
    bt_dict = {"lower": bt_lower, "upper": bt_upper, "middle": bt_middle}
    highest_score, start_key = find_score_and_matrix_to_start_backtrack(
        score_lower, score_upper, score_middle
    )
    first_align, second_align = do_backtrack(bt_dict, start_key, first, second)
    print(highest_score)
    print("".join(first_align))
    print("".join(second_align))


if __name__ == "__main__":
    Fire(main)
