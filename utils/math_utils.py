"""
Getting down to the raw arithmetic
"""

def expected_number_of_kmers_over_multiple_sequences(k, L, N):
    """

    :param k: size of the k-mer nucleotide
    :param L: length of each random DNA string
    :param N: number of random DNA strings
    :return:
    """

    return ( (0.25 ** k) * (L - k + 1) ) * N