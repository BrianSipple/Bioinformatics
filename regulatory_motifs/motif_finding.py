"""
A regulatory motif is a pattern that appears at least once (perhaps with variation)
in each of many different regions that are scattered throughout the genome.

As such, we're now living in a completely different context then when we were
trying to find a DNA box (e.g. finding the Origin or Replication) -- because
a DnaA box is a pattern that clumps, or appears frequently,
within a relatively short interval of the genome.

Concretely, we're now computing over DNA arrays where each item in the
array will represent a `slice` of the genome at a point that's relatively
distant from each other slice.
"""
from math import log
from data_structures.arrays import FrequencyArrays
from origin_of_replication.pattern_finding import neighbors, pattern_count, hamming_distance


def find_motifs(k, d, dna_seqs):
    """
    Finds all (k,d)-motifs in a collection of strings

    A k-mer is a (k,d)-motif if it appears in every
    string from Dna with at most d mismatches

    :param k: int - the length of the k-mer to find
    :param d: int - the number of maximum mistmatches
    :param dna_seqs: a list of DNA string from different positions
    of a genome

    :return: list of strings - all (k, d)-mers
    """
    neighborhood = []
    motifs = []

    ## Generate all k-mer neighbors with at most d-mismatches for each dna_seq
    for i in range(len(dna_seqs)):

        sequence_neighbors = []

        for j in range(len(dna_seqs[i]) - k):

            kmer = dna_seqs[i][j:j+k]
            kmer_neighbors = neighbors(kmer, d)

            for l in range(len(kmer_neighbors)):
                if kmer_neighbors[l] not in sequence_neighbors:
                    sequence_neighbors.append(kmer_neighbors[l])

        neighborhood.append(sequence_neighbors)


    for i in range(len(neighborhood)):
        neighbor_list = neighborhood[i]  # each neighbor list will correspond to one of the passed-in DNA seqs

        for n in neighbor_list:

            is_motif = True

            for seq in range(len(dna_seqs)):
                if pattern_count(dna_seqs[seq], n, mutation_thresh=d) > 0:
                    continue
                else:
                    is_motif = False

            if is_motif:
                if n not in motifs:
                    motifs.append(n)

    return motifs



def motif_matrix_score(motifs):
    """
    When we've constructed a t x k matrix of k-mer motifs
    in t sequences, we can compute the score as the number of
    unpopular letters (nucleotides not matching the most popular nucleotide)
    in each column

    :param motifs: list of lists - k x t k-mer sequences in a motif matrix
    :return: number of unpopular nucleotides in the matrix
    """
    pass




def profile_matrix(motifs):
    """
    Computes frequencies for each nucleotide at each column in `motif`

    :param motifs: list of lists - k x t k-mer sequences in a motif matrix
    :return: 4-k matrix of frequencies for A, C, G, and T in each column
    """
    pass



def consensus_string(motifs):
    """
    Assembles a string of nucleotides based upon the highest-frequency
    nucleotide for each column of `motifs`

    :param motifs: list of lists - k x t k-mer sequences in a motif matrix
    :return: String of nucleotides representing an ideal candidate regulatory motif
    """

def matrix_entropy(matrix):
    """

    :param matrix:
    :return:
    """
    column_entropies = []

    profile = profile_matrix(matrix)

    freq_lists = list(profile.values())

    for i in range(len(freq_lists)):

        column_entropy = 0

        for l in freq_lists:
            column_entropy += (-1 * (l[i] * log(l[i], base=2)))

        column_entropies.append(column_entropy)

    return column_entropies



def median_string(k, dna_strings):
    """
    Finds a k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern
    :param k:
    :param dna_strings:
    :return:
    """
    current_distance = len(dna_strings)**len(dna_strings)  # Set a higher-than-possible max distance
    median = None

    for i in range(4**k - 1):
        pattern = FrequencyArrays._number_to_pattern(i, k)

        m_distance = matrix_distance(pattern, dna_strings)

        if current_distance > m_distance:
            current_distance = m_distance
            median = pattern

    return median


def matrix_distance(pattern, dna_strings):
    """
    Computes the total of the minimum hamming distances of a pattern and dna_string
    over all strings in `dna_strings`
    """

    k = len(pattern)
    distance = 0

    for dna_string in dna_strings:

        found_hamming_distance = len(dna_string) ** len(dna_strings)  # Initialize a maximum

        for i in range(len(dna_string) - k):

            dna_kmer = dna_string[i: i + k]
            hd = hamming_distance(dna_kmer, pattern)

            if found_hamming_distance > hd:
                found_hamming_distance = hd

        distance += found_hamming_distance

    return distance





