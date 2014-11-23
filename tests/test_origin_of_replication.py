import logging
import unittest
from origin_of_replication.pattern_finding import find_clumps, \
    compute_gc_skew, pattern_count, most_frequent_kmers, find_pattern_positions, reverse_complement, hamming_distance, \
    neighbors


class OriginOfReplicationTest(unittest.TestCase):
    def setUp(self):
        super(OriginOfReplicationTest, self).setUp()

    def tearDown(self):
        super(OriginOfReplicationTest, self).tearDown()


    def test_pattern_count(self):
        DNA = "GCGCG"
        pattern = "GCG"

        with self.assertRaises(ValueError):
            pattern_count(DNA, pattern, -1)

        with self.assertRaises(ValueError):
            pattern_count(DNA, pattern, 10)

        self.assertEqual(pattern_count(DNA, pattern), 2)

        DNA = "AAAAAA"
        pattern = "A"

        self.assertEqual(pattern_count(DNA, pattern), 6)

        DNA = "AGCCAGATGGGACTTAGAT"
        pattern = "AGAT"

        self.assertEqual(pattern_count(DNA, pattern), 2)

        DNA = "AAGAGGTTACCCATGTC"
        pattern = "GGG"

        self.assertEqual(pattern_count(DNA, pattern), 0)



    def test_most_frequent_kmers(self):

        DNA = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        K = 4

        self.assertEqual(most_frequent_kmers(DNA, K, 0), set(['CATG', 'GCAT']))



    def test_most_frequent_kmers_with_mismatches(self):

        DNA = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        K = 4
        mutation_thresh = 1
        kmers = most_frequent_kmers(DNA, K, mutation_thresh=mutation_thresh)
        self.assertEqual(kmers, set(["GATG", "ATGC", "ATGT"]))

        # DNA = "CAAAGCAGCTTTCTCTGGGAGCTGGGTTTTGGGTGGGAGCTTTCAACAACTCCTCAGCCTCTTTAGCTTTCTCTGGGCAATGGGCTCAGCCTCCTCCAAAGCTTTTTTCTCCAACAATTTAGCAGCCAAAGCTGGGTGGGCAATTTTGGGAGCCAACTCTGGGTTTTTTAGCTTTCTCCAACAATGGGCTCCTCAGCAGCAGCAGCTGGGAGCAGCTTTTGGGCAACTCAGCCAATTTAGCCTCCTCTGGGCAATGGGAGCTGGGTGGGCAACTCCAATGGGTTTAGCAGCTGGGCAACAACAAAGCTGGGTGGGCTCTGGGCTCTTTAGCAGCAGCCAATTTTGGGCTC"
        # K = 8
        # mutation_thresh = 2
        # kmers = most_frequent_kmers(DNA, K, mutation_thresh)
        # self.assertEqual(kmers, set(["TCAGCAAC"]))



    def test_most_frequent_kmers_with_mismatches_and_reverse_complements(self):

        DNA = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
        K = 4
        mutation_thresh = 1
        kmers = most_frequent_kmers(DNA, K, mutation_thresh, reverse=True)
        self.assertEqual(kmers, set(["ATGT", "ACAT"]))

        DNA = "ATTTGATCATATCATCATTTACATATATTTGTTGATCATTTGATCATTATTATTACTACATCATTACATTACTTGATTTTGATCATTTACATTTGTTGTACATCATATCATATCATTTGATATCATTATATCATTATATTTTGATATTATCTACATATTACATATATCATTACTACATCATTTTGATTATCATTTGTTGAT"
        K = 10
        mutation_thresh = 2
        kmers = most_frequent_kmers(DNA, K, mutation_thresh=mutation_thresh, reverse=True)
        self.assertEqual(kmers, set(['ATTATATCAT', 'ATAATATGAT', 'ATCATATTAT', 'ATGATATAAT']))



    def test_hamming_distance(self):

        p = "AAACCCTTTGGG"
        q = "AAACCCTTTGGG"
        hd = hamming_distance(p, q)
        self.assertEqual(0, hd)


        p = "AAACCCTTTGGG"
        q = "AAACCCTTTGGC"
        hd = hamming_distance(p, q)
        self.assertEqual(1, hd)


        p = "GGGCCGTTGGT"
        q = "GGACCGTTGAC"
        hd = hamming_distance(p, q)
        self.assertEqual(3, hd)


        p = "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
        q = "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
        hd = hamming_distance(p, q)
        self.assertEqual(23, hd)


    def test_find_pattern_count_allowing_for_mutations(self):

        DNA = "CGTGACAGTGTATGGGCATCTTT"
        pattern = "TGT"
        m = 1
        count = pattern_count(DNA, pattern, start=0, end=0, mutation_thresh=m)
        self.assertEqual(8, count)

        DNA = "AACAAGCTGATAAACATTTAAAGAG"
        pattern = "AAAAA"
        m = 2
        count = pattern_count(pattern=pattern, DNA=DNA, mutation_thresh=m)
        self.assertEqual(11, count)


        DNA = "TTTAGAGCCTTCAGAGG"
        pattern = "GAGG"
        m = 2
        count = pattern_count(DNA=DNA, pattern=pattern, mutation_thresh=m)
        self.assertEqual(count, 4)


        DNA = "TTACCGGGCATAGTACTCATGGCTGCCGTTTGAGCCCTCGCTTACTCTATCTTGTAATGATACCAGCACAATGCTTTCTGTCCAAGACCTGATCGGACTATTAGACCGCGACCAGTATTGTGATACCCGTTAGCGCGAAACGACGCACGAATGGGGAAATCCGCCGGATGAAACTGTAAAATTGGCGGGCGCACTCAATTTGTGGACGGACTTCTTGGTAGTTTCTCCTTAAGTCAGGGCTGGCATGAATGAATAGCGTCATGGGCCGATTTACATCAGAATCACAACTGAGTCAGACCAGGCAG"
        pattern = "AATTTGT"
        m = 2
        count = pattern_count(DNA=DNA, pattern=pattern, mutation_thresh=m)
        self.assertEqual(count, 6)



    def test_find_positiions_of_patterns_allowing_for_mismatches(self):

        pattern = "ATTCTGGA"
        DNA = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
        mutation_thresh = 3

        positions = find_pattern_positions(pattern, DNA, mutation_thresh)
        self.assertEqual(positions, [6, 7, 26, 27])



    def test_find_clumps(self):

        # We'll warmup a small test first
        DNA = "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
        k = 5
        L = 50
        t = 4
        res = {'CGACA', 'GAAGA'}
        self.assertEqual(find_clumps(DNA, k, L, t), res)

        DNA = "CCGACAGGCTAGTCTATAATCCTGAGGCGTTACCCCAATACCGTTTACCGTGGGATTTGCTACTACAACTCCTGAGCGCTACATGTACGAAACCATGTTATGTAT"
        k = 4
        L = 30
        t = 3
        self.assertEqual(find_clumps(DNA, k, L, t), set(['CTAC', 'TACC', 'ATGT']))

        # For when we want to test on ecoli...
        # k = 9
        # L = 500
        # t = 3
        # print("Clump test...")
        #
        # with open("../resources/ecoli.txt", 'r') as file:
        #     ecoli = file.read()
        #
        # clumps = find_clumps(ecoli, k, L, t)
        #
        # self.assertEquals(clumps, ("CGACA", "GAAGA"))



    def test_compute_gc_skew(self):

        genome = "GAGCCACCGCGATA"
        true_skew = [0, 1, 1, 2, 1, 0, 0, -1, -2, -1, -2, -1, -1, -1, -1]
        true_min_skew_positions = [8, 10]

        test_skew, test_min_skew_positions = compute_gc_skew(genome, chart=False)

        self.assertEqual(true_skew, test_skew)
        self.assertEqual(true_min_skew_positions, test_min_skew_positions)

        # genome = "GATACACTTCCCAGTAGGTACTG"
        # skew_array, min_skew_positions = compute_gc_skew(genome)
        #
        # print(skew_array, min_skew_positions)


    def test_reverse_complement(self):

        DNA = "TTGTGTC"
        r_comp = reverse_complement(DNA, as_string=True)
        self.assertEqual(r_comp, "GACACAA")



    def test_neighbors(self):

        pattern = "ACGT"
        d = 3

        n = neighbors(pattern, d)
        self.assertEqual(len(n), 175)
        print(n)



if __name__ == "__main__":
    unittest.main()

