import unittest
from math_utils import expected_number_of_kmers_over_multiple_sequences
from regulatory_motifs.motif_finding import find_motifs, matrix_distance, median_string


class RegulatoryMotifsTest(unittest.TestCase):
    def setUp(self):
        super(RegulatoryMotifsTest, self).setUp()

    def tearDown(self):
        super(RegulatoryMotifsTest, self).tearDown()


    def test_expected_value_over_multiple_sequences(self):
        k = 9
        L = 1000
        N = 500

        expected_number_of_kmers = expected_number_of_kmers_over_multiple_sequences(k, L, N)

        self.assertEqual(1.89208984375, expected_number_of_kmers)


    def test_finding_motifs(self):
        k = 3
        d = 1
        dna_seqs = [
            "ATTTGGC",
            "TGCCTTA",
            "CGGTATC",
            "GAAAATT"
        ]
        #
        # with open('...', 'r') as f:
        # dna_seqs.append(f.readline())

        motifs = find_motifs(k, d, dna_seqs)
        self.assertListEqual(["ATA", "ATT", "GTT", "TTT"], motifs)

        #print(" ".join([m for m in motifs]))

    def test_matrix_distance(self):
        pattern = "AAA"
        dna_strings = ["TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"]

        m_dist = matrix_distance(pattern, dna_strings)
        self.assertEqual(5, m_dist)


    def test_median_string(self):
        k = 3
        dna_strings = [
            "AAATTGACGCAT",
            "GACGACCACGTT",
            "CGTCAGCGCCTG",
            "GCTGAGCACCGG",
            "AGTACGGGACAG"
        ]

        self.assertEqual(median_string(k, dna_strings), "GAC")


if __name__ == "__main__":
    unittest.main()