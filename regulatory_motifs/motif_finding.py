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
from origin_of_replication.pattern_finding import neighbors, pattern_count


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



