import unittest
from utils import hammingDistance

class OriginOfReplicationTest(unittest.TestCase):

	def setUp(self):
		super(OriginOfReplicationTest, self).setUp()

	def tearDown(self):
		super(OriginOfReplicationTest, self).tearDown()

	
	def test_hamming_distance(self):

		p = "AAACCCTTTGGG"
		q = "AAACCCTTTGGC"

		hd = hammingDistance(p, q)

		self.assertEquals(1, hd)



if __name__ == "__main__":

	unittest.main()

