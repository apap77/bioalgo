import sys; sys.path.append('../')
import bioalgo.genome as genome
import unittest

class AssemblyTest(unittest.TestCase):

	def setUp(self):
		self.db = genome.assembly.DeBrujin('TAATGCCATGGGATGTT', k=3)

	def test_de_brujin(self):
		
		self.assertEqual(self.db['TA'], ['AA'])
		self.assertEqual(self.db['AT'], ['TG', 'TG', 'TG'])
		self.assertEqual(self.db['TG'], ['GC', 'GG', 'GT'])

	def test_de_brujin_str(self):
		s = 'AA: AT\nAT: TG, TG, TG\nCA: AT\nCC: CA\nGA: AT\nGC: CC\nGG: GG, GA\nGT: TT\nTA: AA\nTG: GC, GG, GT'
		self.assertEqual(str(self.db), s)