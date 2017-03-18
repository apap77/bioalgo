import sys; sys.path.append('../')
import bioalgo.align as align
import unittest

class AlignTest(unittest.TestCase):

	def setUp(self):
		self.aa1 = 'PLEASANTLY'
		self.aa2 = 'MEANLY'
		self.dna1 = 'AACCTTGG'
		self.dna2 = 'ACACTGTGA'

	def test_basic_global_aligner_assertion(self):
		bga = align.BasicGlobalAligner(self.aa1, self.aa2)

		self.assertRaises(ValueError, print, bga)

	def test_basic_global_alinger_aa_align(self):
		bga = align.BasicGlobalAligner(self.aa1, self.aa2)
		bga.align()

		self.assertEqual(str(bga), 'PLEASANTLY\n-MEA--N-LY')

	def test_basic_global_alinger_longest_multiple_subsequence(self):
		bga = align.BasicGlobalAligner(self.dna1, self.dna2, scoringMatrix='nuc_identity', indelPenalty=0)
		bga.align()

		self.assertEqual(str(bga), 'A-ACCT-TG-G\nACAC-TGTGA-')

class ScoringMatrixTest(unittest.TestCase):

	def test_blosum62(self):
		self.assertEqual(align.ScoringMatrix('blosum62')['A']['A'], 4)

	def test_nuc_blast(self):
		self.assertEqual(align.ScoringMatrix('nuc_blast')['A']['A'], 5)

	def test_nuc_identity(self):
		self.assertEqual(align.ScoringMatrix('nuc_identity')['A']['A'], 1)
		self.assertEqual(align.ScoringMatrix('nuc_identity')['A']['C'], 0)

	def test_pam250(self):
		self.assertEqual(align.ScoringMatrix('pam250')['A']['A'], 2)