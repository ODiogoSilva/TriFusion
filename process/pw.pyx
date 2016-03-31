
def calc_similarity(seq1, seq2):
    """
    Gets the similarity between two sequences in proportion. The
    proportion is calculated here so that sites with only missing
    data/gaps are not accounted in the total length.
    :param seq1: string
    :param seq2: string
    """
    similarity = 0.0
    effective_len = 0.0
    missing = ["n", "-"]
    for c1, c2 in zip(*[seq1, seq2]):
        # Ignore comparisons with ONLY missing data / gaps
        if c1 in missing or c2 in missing:
            continue
        elif c1 == c2:
            similarity += 1.0
            effective_len += 1.0
        else:
            effective_len += 1.0
    if effective_len:
        return similarity, effective_len
    else:
        return None, None