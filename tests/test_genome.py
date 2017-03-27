import sys; sys.path.append('../')
import bioalgo.genome as genome
import unittest
class AssemblyTest(unittest.TestCase):

	def setUp(self):
		patterns = ['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG']
		patterns2 = ['CTTA', 'ACCA', 'TACC', 'GGCT', 'GCTT', 'TTAC']
		self.dbText = genome.assembly.DeBrujin('TAATGCCATGGGATGTT', k=3)
		self.dbText2 = genome.assembly.DeBrujin('TAATGCCATGGGATGTT', k=4)
		self.dbPattern = genome.assembly.DeBrujin(patterns, k=4)
		self.dbPattern2 = genome.assembly.DeBrujin(patterns2, k=4)

	def test_de_brujin(self):
		self.assertEqual(self.dbText['TA'], ['AA'])
		self.assertEqual(self.dbText['AT'], ['TG', 'TG', 'TG'])
		self.assertEqual(self.dbText['TG'], ['GC', 'GG', 'GT'])

		self.assertEqual(self.dbPattern['AGG'], ['GGG'])
		self.assertEqual(self.dbPattern['GGG'], ['GGG', 'GGA'])

	def test_de_brujin_str(self):
		s = 'AA -> AT\nAT -> TG, TG, TG\nCA -> AT\nCC -> CA\nGA -> AT\nGC -> CC\nGG -> GG, GA\nGT -> TT\nTA -> AA\nTG -> GC, GG, GT'
		self.assertEqual(str(self.dbText), s)

	def test_de_brujin_reconstruct(self):
		self.assertEqual(self.dbText.reconstruct(), 'TAATGCCATGGGATGTT')
		self.assertEqual(self.dbText2.reconstruct(), 'TAATGCCATGGGATGTT')
		self.assertRaises(AssertionError, self.dbPattern.reconstruct)
		self.assertEqual(self.dbPattern2.reconstruct(), 'GGCTTACCA')