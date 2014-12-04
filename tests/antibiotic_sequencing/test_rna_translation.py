import unittest
from antibiotic_sequencing.rna_translation import translate_to_peptides, compute_possible_rna_origins, \
    compute_possible_dna_origins, compute_cyclopeptides, compute_mass_spectrum, count_peptides_with_mass, \
    compute_peptide_total_mass, count_subpeptides_in_linear_peptide, count_subpeptides_in_cyclopeptides, \
    find_cyclopeptide_in_mass_spectrum
from utils import spectrum_to_peptides


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
        with open('../input_data/bacillus_brevis.txt', 'r') as f:
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


    def test_peptide_mass(self):

        peptide = "NQEL"
        total_mass = compute_peptide_total_mass(peptide)

        self.assertEqual(484, total_mass)



    def test_count_subpeptides_in_linear_peptide(self):

        peptide_length = 5
        self.assertEqual(17, count_subpeptides_in_linear_peptide(peptide_length))

        peptide_length = 23250
        #print(count_subpeptides_in_linear_peptide(peptide_length))


    def test_count_subpeptides_in_cyclopeptide(self):

        peptide_length = 4
        expected = peptide_length * (peptide_length-1) + 2
        self.assertEqual(expected, count_subpeptides_in_cyclopeptides(peptide_length))


    def test_find_cyclopeptide_in_experimental_spectrum(self):

        experimental_spectrum = "0 71 101 113 131 184 202 214 232 285 303 315 345 416"
        experimental_spectrum = list(map(int, experimental_spectrum.split(" ")))  # convert to array of ints

        results = find_cyclopeptide_in_mass_spectrum(experimental_spectrum)

        output = ""

        for res in results:
            output += ("-".join(map(str, res))) + " "

        print(output)

        for result in results:
            print(spectrum_to_peptides(result))


    def pending_test_convolution_cyclopeptide_sequencing(self):

        m = 20
        n = 60
        exp_spec = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]

        print(convolution_cyclopeptide_sequencing(m, n, exp_spec))


if __name__ == "__main__":
    unittest.main()