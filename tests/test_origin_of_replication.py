import logging
import unittest
from origin_of_replication.pattern_finding import neighbors, find_approximate_pattern_positions, find_clumps, \
    compute_gc_skew, pattern_count, most_frequent_kmers
from utils import hammingDistance, make_spaced_string_from_comma_separated_array


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

        self.assertEqual(most_frequent_kmers(DNA, K), set(['CATG', 'GCAT']))



    def test_hamming_distance(self):
        p = "AAACCCTTTGGG"
        q = "AAACCCTTTGGG"

        hd = hammingDistance(p, q)
        self.assertEqual(0, hd)

        p = "AAACCCTTTGGG"
        q = "AAACCCTTTGGC"

        hd = hammingDistance(p, q)
        self.assertEqual(1, hd)

        p = "GGGCCGTTGGT"
        q = "GGACCGTTGAC"

        hd = hammingDistance(p, q)
        self.assertEqual(3, hd)


    def test_neighbors(self):
        pattern = "ACG"
        d = 1

        neighborhood = neighbors(pattern, d)

        self.assertEqual(neighborhood,
            [
                "CCG",
                "TCG",
                "GCG",
                "AAG",
                "ATG",
                "AGG",
                "ACA",
                "ACC",
                "ACT",
                "ACG"
            ]
        )

        # pattern = "GCAAAAGCGT"
        # d = 3
        # print(make_spaced_string_from_comma_separated_array(neighbors(pattern, d)))


    def _test_find_apporoximate_pattern_matches(self):
        pattern = "ATTCTGGA"
        genome = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
        d = 3

        positions = find_approximate_pattern_positions(pattern, genome, d)

        self.assertEqual(positions, [6, 7, 26, 27])


    def test_find_clumps(self):
        k = 9
        L = 500
        t = 3

        with open("../resources/ecoli.txt", 'r') as file:
            ecoli = file.read()

        clumps = find_clumps(ecoli, k, L, t)

        self.assertEquals(clumps, ("CGACA", "GAAGA"))


    def test_compute_gc_skew(self):
        genome = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
        skew = compute_gc_skew(genome, chart=True)


if __name__ == "__main__":
    unittest.main()

