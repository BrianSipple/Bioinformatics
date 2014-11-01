import unittest
from antibiotic_sequencing.rna_translation import find_rna_translation


class RNATranslationTest(unittest.TestCase):
    def setUp(self):
        super(RNATranslationTest, self).setUp()

    def tearDown(self):
        super(RNATranslationTest, self).tearDown()


    def test_rna_translation(self):

        RNA = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
        peptides = find_rna_translation(RNA)

        print(peptides)
        self.assertEqual("MAMAPRTEINSTRING", peptides)

        RNA = ""



if __name__ == "__main__":
    unittest.main()