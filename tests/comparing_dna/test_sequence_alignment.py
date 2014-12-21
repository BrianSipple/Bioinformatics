import unittest
from comparing_dna.sequence_alignment import min_number_of_coins


class SequenceAlignmentTest(unittest.TestCase):

    def setUp(self):
        super(SequenceAlignmentTest, self).setUp()

    def tearDown(self):
        super(SequenceAlignmentTest, self).tearDown()


    def test_min_number_of_coins(self):

        money = 19914
        coin_choices = [20,19,15,13,9,5,3,1]

        min_num = min_number_of_coins(money, coin_choices)

        self.assertEqual(min_num, 996)


if __name__ == "__main__":
    unittest.main()