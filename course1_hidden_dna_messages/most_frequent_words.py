from collections import Counter

from fire import Fire

def most_frequent_words(text, k):
    max_freq = 0
    max_freq_kmer = []
    kmer_map = {}
    for ix in range(len(text)-k+1):
        cur_kmer = text[ix:ix+k]
        if cur_kmer in kmer_map:
            kmer_map[cur_kmer] += 1
        else:
            kmer_map[cur_kmer] = 1
        cur_cnt = kmer_map[cur_kmer]
        if cur_cnt > max_freq:
            max_freq_kmer = [cur_kmer]
            max_freq = cur_cnt
        elif cur_cnt == max_freq:
            max_freq_kmer.append(cur_kmer)
    return max_freq_kmer


if __name__ == "__main__":
    Fire(most_frequent_words)
