import unittest
from antibiotic_sequencing.rna_translation import translate_to_peptides, compute_possible_rna_origins, \
    compute_possible_dna_origins


class RNATranslationTest(unittest.TestCase):
    def setUp(self):
        super(RNATranslationTest, self).setUp()

    def tearDown(self):
        super(RNATranslationTest, self).tearDown()



    def test_rna_translation(self):

        RNA = "AUG"
        peptide = translate_to_peptides(RNA)
        self.assertEqual("M", peptide)

        RNA = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        peptides = translate_to_peptides(RNA)
        self.assertEqual("MAMAPRTEINSTRING", peptides)



    def test_compute_possible_dna_origins(self):

        peptides = ['V', 'K', 'L', 'F', 'W', 'P', 'F', 'N', 'Q', 'Y']
        count = compute_possible_rna_origins(peptides)
        self.assertEqual(count, 27)



    def test_compute_possible_dna_origins(self):

        DNA = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
        peptides = "MA"

        possible_origins = compute_possible_dna_origins(DNA, peptides)
        print(possible_origins)



if __name__ == "__main__":
    unittest.main()