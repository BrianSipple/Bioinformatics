"""
Optimized Data structures based on arrays
"""
import itertools


class FrequencyArrays(object):
    """
    Construct a frequency array for a given array string
    """
    def init(self, DNA, K):
        """
        Constructs a frequency array for the given DNA and K

        :param DNA: String - DNA string
        :param K: Integer - Length of the K-mer
        """
        assert K >= 1
        assert len(DNA) >= K

        self.K = K

        # All possible combinations of k-mers
        self.kmers = [''.join(p) for p in itertools.product(['A', 'C', 'T', 'G'], repeat=K)]

        # Prepare frequency array
        self.freq = [0] * pow(4, K)

        #####################
        # Precompute dictionaries of pattern_to_number and number_to_pattern
        ####################
        self.pattern_to_number = {
            kmer: self._pattern_to_number(kmer) for kmer in self.kmers
        }
        self.number_to_pattern = {
            i: self._number_to_pattern(i, K) for i in range(len(self.freq))
        }

        # compute frequency array
        for i in range(len(DNA) - K + 1):
            index = self.pattern_to_number[''.join(DNA[i : i + K])]
            self.freq[index] += 1


    def get_frequency_array(self):
        """
        Return already computed frequency array
        """
        return self.freq


    def get_most_frequent(self):
        """
        Return a list with the most frequent k-mers
        :return:
        """
        max_count = max(self.freq)
        most_freq = set()

        for i in range(len(self.kmers)):
            if self.freq[i] == max_count:
                most_freq.add(self.number_to_pattern[i, self.K])
        return most_freq


    def set_frequency(self, kmer, f):
        """
        Set f as the frequency of the kmer in the frequency array.

        This method allows the use of frequency array for sliding windows in
         a DNA string
        :param kmer: String - kmer for which to set frequency
        :param f: Integer - frequency
        """
        assert f >= 0
        self.freq[self.pattern_to_number[kmer]] = f


    def _pattern_to_number(self, kmer):
        """
        Returns the index of kmer in the frequency array.

        The approach to computing _pattern_to_number(kmer) is based on a simple
        observation. If we remove the final symbol from all lexicographically
        ordered k-mers, the resulting list is still ordered lexicographically.
        In the case of DNA strings, every (k - 1)-mer in the resulting list is
        repeated four times.

        Thus, the number of 3-mers occurring before AGT is equal to four times
        the number of 2-mers occurring before AG plus the number of 1-mers occuring
        before T. Therefore,

        _pattern_to_number(AGT) = 4 * _pattern_to_number(AG) + _symbol_to_number(T)
                                = 8 + 3 = 11
        """
        if len(kmer) == 1:
            return self._symbol_to_number(kmer)
        else:
            return 4 * self._pattern_to_number(kmer[:-1]) + self._symbol_to_number(kmer[-1])


    def _number_to_pattern(self, n, k):
        """
        Returns the k-mer that is represented by n

        In order to compute the inverse function NumberToPattern(index, k), we return to the same equation.

        The equation implies that when we divide n = _pattern_to_number(kmer) by 4,
        the remainder will be equal to _symbol_to_number(s), and the quotient will
        be equal to pattern_to_number(kmer[:]). Thus, we can use this fact to peel
        away symbols at the end of kmer one at a time.
        """
        if k == 1:
            return self._number_to_symbol(k)
        else:
            quotient, remainder = divmod(n, 4)
            return self._number_to_pattern(quotient, k-1) + self._number_to_symbol(remainder)


    def _symbol_to_number(self, symbol):
        """
        Function transforming symbols A, C, G, and T into the respective integers 0, 1, 2, and 3
        """
        nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        return nucleotides[symbol]


    def _number_to_symbol(self, n):
        """
        Function transforming symbols o, 1, 2, and 3 into the respective nucleotides A, C, G, and T.
        """
        nucleotides = ['A', 'C', 'G', 'T']
        return nucleotides[n]




