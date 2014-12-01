
"""
Getting down to the raw arithmetic
"""
from math import log


def expected_number_of_kmers_over_multiple_sequences(k, L, N):
    """

    :param k: size of the k-mer nucleotide
    :param L: length of each random DNA string
    :param N: number of random DNA strings
    :return:
    """

    return ( (0.25 ** k) * (L - k + 1) ) * N


def entropy(distribution):
    """
    Calculates entropy for a probability distribution
    """
    return -1 * sum([ (prob * log(prob, 2)) if prob != 0 else 0 for prob in distribution ])