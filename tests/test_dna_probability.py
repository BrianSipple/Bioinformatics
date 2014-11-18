import unittest
from molecular_clocks.dna_probability import expected_number_of_kmers_over_multiple_sequences


class DNAProbabilityTest(unittest.TestCase):

    def setUp(self):
        super(DNAProbabilityTest, self).setUp()

    def tearDown(self):
        super(DNAProbabilityTest, self).tearDown()



    def test_expected_value_over_multiple_sequences(self):

        k = 9
        L = 1000
        N = 500

        expected_number_of_kmers = expected_number_of_kmers_over_multiple_sequences(k, L, N)

        self.assertEqual(1.89208984375, expected_number_of_kmers)



if __name__ == "__main__":
    unittest.main()