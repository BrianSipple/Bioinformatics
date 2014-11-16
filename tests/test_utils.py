import unittest
from utils import spectrum_to_peptides


class UtilsTest(unittest.TestCase):
    def setUp(self):
        super(UtilsTest, self).setUp()

    def tearDown(self):
        super(UtilsTest, self).tearDown()


    def test_spectrum_to_peptides(self):

        unique_peptides = ['K', 'P', 'V', 'M', 'T', 'A', 'E', 'H', 'C', 'I', 'N', 'D', 'F', 'S', 'G', 'Y', 'R', 'W']
        unique_masses = [128, 97, 99, 131, 101, 71, 129, 137, 103, 113, 114, 115, 147, 87, 57, 163, 156, 186]

        self.assertEqual(unique_peptides, spectrum_to_peptides(unique_masses))

        IND_masses = [113, 114, 115]
        self.assertEqual(["I", "N", "D"], spectrum_to_peptides(IND_masses))





if __name__ == "__main__":

    unittest.main()

