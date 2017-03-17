import sys; sys.path.append('../')
import bioalgo.seq as seq
import unittest

class SeqTest(unittest.TestCase):

	def test_dna_nucleotide_assertion(self):
		self.assertRaises(AssertionError, seq.DNA, 'ATGCU')

	def test_rna_nucleotide_assertion(self):
		self.assertRaises(AssertionError, seq.RNA, 'ATGCU')

if __name__ == '__main__':
	d = seq.RNA('AUGAAUGAUUGA')
	kc = d.kmer_composition(k=3)
	print(kc)
	print(d.translate(frameStart=-1))
	print(d.reverse_complement())

	d = seq.DNA('ATGTTATGGATAG')
	kc = d.kmer_composition(k=2)
	print(kc)