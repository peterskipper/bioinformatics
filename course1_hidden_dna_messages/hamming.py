from itertools import islice

from fire import Fire


def hamming_dist(str1, str2):
    dist = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            dist += 1
    dist += abs(len(str1) - len(str2))
    return dist


def window(seq, n=2):
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def count_2(text, substr):
    cnt = 0
    n = len(substr)
    for wnd in window(text, n):
        if hamming_dist(wnd, substr) <= 2:
            cnt += 1
    return cnt


if __name__ == "__main__":
    Fire(count_2)
