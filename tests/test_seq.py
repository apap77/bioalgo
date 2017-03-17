import sys; sys.path.append('../')
import bioalgo.seq as seq
from bioalgo.seq import codon
import unittest

class SeqTest(unittest.TestCase):

	def test_dna_nucleotide_assertion(self):
		self.assertRaises(AssertionError, seq.DNA, 'ATGCU')

	def test_rna_nucleotide_assertion(self):
		self.assertRaises(AssertionError, seq.RNA, 'ATGCU')

	def test_reverse_complement(self):
		d = seq.DNA('ATGC')
		r = seq.RNA('AUGC')

		self.assertEqual(d.reverse_complement(), 'GCAT')
		self.assertEqual(r.reverse_complement(), 'GCAU')

	def test_transcribe(self):
		d = seq.DNA('ATGC')
		self.assertEqual(d.transcribe(), 'GCAU')
		self.assertEqual(d.transcribe(reverse=False), 'UACG')

	def test_translate(self):
		r = seq.RNA('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
		self.assertEqual(r.translate(), 'MAMAPRTEINSTRING*')

	def test_kmer_composition(self):
		d = seq.DNA('ATTGCCATTGGC')
		self.assertEqual(d.kmer_composition(k=3)['ATT'], 2)

class CodonTest(unittest.TestCase):

	def setUp(self):
		self.c = codon.CodonTable()

	def test_codon(self):
		self.assertEqual(self.c['UUU'], 'F')

	def test_stop_codon(self):
		self.assertEqual(self.c['UGA'], '*')
		self.assertEqual(self.c['UAG'], '*')
		self.assertEqual(self.c['UAA'], '*')

	def test_get_codons(self):
		self.assertEqual(set(self.c.get_codons('T')), set(['ACU', 'ACC', 'ACA', 'ACG']))