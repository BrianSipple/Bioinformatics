"""
Data structures based on dictionaries
"""

from collections import defaultdict


class FrequencyDict(defaultdict):
    """
    Data structure to store frequency of k-mer in a DNA string
    """

    def __init__(self, DNA, K, m=0):
        super(FrequencyDict, self).__init__(int)

        # Increment the frequency of all integers
        for i in range(len(DNA) - K + 1):
            kmer = DNA[i: i + K]
            self[kmer] += 1

            if m:
                mutants = set()
                self._mutate(mutants, kmer, kmer, K - 1, m)
                mutants.remove(kmer)
                for kmer in mutants:
                    self[kmer] += 1


    def _hamming(self, s1, s2):
        """
        Computes the Hamming Distance between s1 and s2
        """
        assert len(s1) == len(s2)
        return sum([1 if c1 != c2 else 0 for c1, c2 in zip(s1, s2)])

    def _mutate(self, mutations, kmer, orig_kmer, pos, m):
        """
        Return a set containing all mutations of orig_kmer with at most m
        mismatches from it.
        """
        if pos >= 0:
            for base in ['A', 'C', 'T', 'G']:
                kmer_mut = list(kmer)
                kmer_mut[pos] = base
                kmer_mut = ''.join(kmer_mut)
                if self._hamming(kmer_mut, orig_kmer) <= m:
                    mutations.add(kmer_mut)
                    self._mutate(mutations, kmer_mut, orig_kmer, pos - 1, m)

