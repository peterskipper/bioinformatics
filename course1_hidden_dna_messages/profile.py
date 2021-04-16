def make_pseudocnt(nucs):
    assert set(nucs).issubset(
        {"A", "C", "G", "T"}
    ), "Expected A/C/G/T nucleotides, found {}".format(set(nucs))
    c = Counter(nucs)
    # make pseudo counts
    for nuc in ["A", "C", "G", "T"]:
        c[nuc] += 1
    return c


def make_percent(cntr):
    ttl = sum(cntr.values())
    return {k: v / ttl for k, v in cntr.items()}


def profile(motifs):
    assert len(set([len(m) for m in motifs])) == 1, "motifs are not all the same length"

    profile = []
    for ix in range(len(motifs[0])):
        nucs = [m[ix] for m in motifs]
        psd_cnts = make_pseudocnt(nucs)
        psd_cnts = make_percent(psd_cnts)
        profile.append(psd_cnts)

    return pd.DataFrame(profile).T
