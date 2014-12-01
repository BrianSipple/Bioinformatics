import unittest
from math_utils import entropy


class MathUtilsTest(unittest.TestCase):

    def setUp(self):
        super(MathUtilsTest, self).setUp()

    def tearDown(self):
        super(MathUtilsTest, self).tearDown()


    def test_entropy(self):

        distribution = [0.5, 0.0, 0.0, 0.5]
        ent = entropy(distribution)

        self.assertEqual(1, ent)

        distribution = [0.0, 0.0, 0.9, 0.1]
        ent = entropy(distribution)

        self.assertEqual(0.47, round(ent, 2))


        distributions = [
            [0.5, 0, 0, 0.5],
            [0.25, 0.25, 0.25, 0.25],
            [0, 0, 0, 1],
            [0.25, 0.0, 0.5, 0.25],
            [0.0, 0.0, 0.9, 0.1],
        ]

        print([entropy(dist) for dist in distributions])


if __name__ == "__main__":
    unittest.main()
