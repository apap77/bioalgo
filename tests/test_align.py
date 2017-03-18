import sys; sys.path.append('../')
import bioalgo.align as align
import unittest

class AlignTest(unittest.TestCase):

	def setUp(self):
		self.aa1 = 'PLEASANTLY'
		self.aa2 = 'MEANLY'

	def test_basic_global_aligner_assertion(self):
		bga = align.BasicGlobalAligner(self.aa1, self.aa2)

		self.assertRaises(ValueError, print, bga)

	def test_basic_global_alinger_align(self):
		bga = align.BasicGlobalAligner(self.aa1, self.aa2)
		bga.align()

		self.assertEqual(str(bga), 'PLEASANTLY\n-MEA--N-LY')