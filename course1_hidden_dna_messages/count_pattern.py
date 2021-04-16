from fire import Fire


def count_pattern(text, pattern):
    count = 0
    pat_size = len(pattern)
    for ix in range(len(text)):
        if text[ix:ix+pat_size] == pattern:
            count += 1
    return f"Count is: {count}"


if __name__ == "__main__":
    Fire(count_pattern)
