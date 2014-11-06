import unittest
from antibiotic_sequencing.rna_translation import translate_to_peptides, compute_possible_rna_origins, \
    compute_possible_dna_origins, compute_cyclopeptides, compute_mass_spectrum


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



    def test_compute_possible_rna_origins(self):

        peptides = ['V', 'K', 'L', 'F', 'W', 'P', 'F', 'N', 'Q', 'Y']
        possible_origins = compute_possible_rna_origins(peptides)
        self.assertEqual(len(possible_origins), 27)



    def test_compute_possible_dna_origins(self):

        DNA = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
        peptides = "MA"
        possible_origins = compute_possible_dna_origins(DNA, peptides)
        self.assertEqual(['ATGGCC', 'ATGGCC', 'GGCCAT'], possible_origins)


        bacillus_brevis = ""
        # bacillus_brevis test
        with open('../resources/bacillus_brevis.txt', 'r') as f:
            while f.readline():
                bacillus_brevis += f.readline().rstrip()

        Tyrocidine_B1 = "VKLFWPFNQY"

        # None exist!
        possible_origins = compute_possible_dna_origins(bacillus_brevis, Tyrocidine_B1)
        self.assertEqual(0, len(possible_origins))



    def test_compute_cyclopeptides(self):

        Tyrocidine_B1 = "VKLFWPFNQY"
        cyclopeptides = compute_cyclopeptides(Tyrocidine_B1)

        #print(cyclopeptides)


    def test_compute_mass_spectrum(self):

        peptide = "LEQN"
        expected_mass_spectrum = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        result = compute_mass_spectrum(peptide)
        self.assertListEqual(expected_mass_spectrum, result)


        peptide = "IAQMLFYCKVATN"
        expected_mass_spectrum = [0, 71, 71, 99, 101, 103, 113, 113, 114, 128, 128, 131, 147, 163, 170, 172, 184, 199, 215, 227, 227, 231, 244, 259, 260, 266, 271, 286, 298, 298, 310, 312, 328, 330, 330, 372, 385, 391, 394, 399, 399, 399, 401, 413, 423, 426, 443, 443, 470, 493, 498, 502, 513, 519, 526, 527, 541, 554, 556, 557, 564, 569, 590, 598, 598, 616, 626, 640, 654, 657, 658, 665, 670, 682, 697, 697, 703, 711, 729, 729, 729, 753, 771, 779, 785, 785, 800, 812, 817, 824, 825, 828, 842, 842, 866, 884, 892, 913, 918, 925, 926, 928, 941, 955, 956, 963, 969, 980, 989, 989, 1012, 1039, 1039, 1056, 1059, 1069, 1081, 1083, 1083, 1088, 1091, 1097, 1110, 1152, 1152, 1152, 1170, 1172, 1184, 1184, 1196, 1211, 1216, 1222, 1223, 1238, 1251, 1255, 1255, 1267, 1283, 1298, 1310, 1312, 1319, 1335, 1351, 1354, 1354, 1368, 1369, 1369, 1379, 1381, 1383, 1411, 1411, 1482]
        result = compute_mass_spectrum(peptide)
        self.assertListEqual(expected_mass_spectrum, result)








if __name__ == "__main__":
    unittest.main()