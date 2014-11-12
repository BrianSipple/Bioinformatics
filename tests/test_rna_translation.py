import unittest
from antibiotic_sequencing.rna_translation import translate_to_peptides, compute_possible_rna_origins, \
    compute_possible_dna_origins, compute_cyclopeptides, compute_mass_spectrum, count_peptides_with_mass, \
    compute_peptide_total_mass, count_linear_subpeptides_in_peptide, count_cyclopeptides_in_peptide


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

        print(cyclopeptides)


    def test_compute_mass_spectrum(self):

        peptide = "NQEL"
        expected_mass_spectrum = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
        result = compute_mass_spectrum(peptide)
        self.assertListEqual(expected_mass_spectrum, result)


        # peptide = "DPHMRVGIHPEHIEQ"
        # result = compute_mass_spectrum(peptide, cyclic=True)
        # print(" ".join([str(num) for num in result]))


    # def test_count_peptides_with_mass(self):
    #
    #     total_mass = 100000
    #     possible_peptides = count_peptides_with_mass(total_mass, recursive=True)
    #
    #     print (possible_peptides)


    def test_peptide_masss(self):

        peptide = "NQEL"
        total_mass = compute_peptide_total_mass(peptide)

        self.assertEqual(484, total_mass)



    def test_count_linear_subpeptides_in_peptide(self):

        peptide_length = 5
        self.assertEqual(120, count_linear_subpeptides_in_peptide(peptide_length))


    def test_count_cyclopeptides_in_peptide(self):

        peptide_length = 5
        expected = peptide_length * (peptide_length-1)

        self.assertEqual(expected, count_cyclopeptides_in_peptide(peptide_length))



if __name__ == "__main__":
    unittest.main()