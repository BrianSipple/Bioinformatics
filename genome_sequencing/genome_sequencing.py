import numpy as np


def sequence_composition(k, sequence):
    """
    Given a string `seqeunce`, its k-mer composition Composition(k, sequence)
    is the collection of all k-mer substrings of sequence (including repeated k-mers)

    :return: the collection of all k-mers in LEXICOGRAPHIC order
    """
    res = []

    for i in range(len(sequence) - k + 1):
        res.append(sequence[i:i+k])

    return list(np.sort(res))