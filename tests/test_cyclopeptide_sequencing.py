import unittest

from antibiotic_sequencing.rna_translation import compute_mass_spectrum, is_linear_spectrum_consistent


class CyclopeptideSequencingTest(unittest.TestCase):
    def setUp(self):
        super(CyclopeptideSequencingTest, self).setUp()

    def tearDown(self):
        super(CyclopeptideSequencingTest, self).tearDown()


    def test_score(self):
        peptide = "PEEP"
        spectrum = [0, 97, 129, 129, 129, 194, 226, 323, 323, 355, 452]
        #print(score(compute_mass_spectrum(peptide, cyclic=False), spectrum))

    def test_consistency(self):

        peptide = "TCE"
        theor_spectrum = [0, 71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231, 298, 303, 328, 330, 332, 333]

        linear_spec = compute_mass_spectrum(peptide, cyclic=False)

        self.assertFalse(is_linear_spectrum_consistent(linear_spec, theor_spectrum))





if __name__ == "__main__":

    unittest.main()

