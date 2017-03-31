import sys; sys.path.append('../')
import bioalgo.seq as seq
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
		self.c = seq.codon.CodonTable()

	def test_codon(self):
		self.assertEqual(self.c['UUU'], 'F')

	def test_stop_codon(self):
		self.assertEqual(self.c['UGA'], '*')
		self.assertEqual(self.c['UAG'], '*')
		self.assertEqual(self.c['UAA'], '*')

	def test_get_codons(self):
		self.assertEqual(set(self.c.get_codons('T')), set(['ACU', 'ACC', 'ACA', 'ACG']))

class ParseTest(unittest.TestCase):

	def test_invalid_fasta(self):
		self.assertRaises(AssertionError, seq.Fasta, './tests/data/invalid.fasta', 'dna')

	def test_dna_fasta(self):
		fastaSequences = seq.Fasta(filePath='./tests/data/dna.fasta', seqType='dna')

		self.assertEqual(len(fastaSequences), 4)
		self.assertEqual(fastaSequences[0].seq[:4], 'GATA')

	def test_invalid_seq_type(self):
		self.assertRaises(AssertionError, seq.Fasta, './tests/data/dna.fasta', 'rna')

	def test_rna_fasta(self):
		fastaSequences = seq.Fasta(filePath='./tests/data/rna.fasta', seqType='rna')

		self.assertEqual(len(fastaSequences), 4)
		self.assertEqual(fastaSequences[0].seq[:4], 'GAUA')

	def test_fastq(self):
		fastqSequences = seq.Fastq(filePath='./tests/data/dna.fastq', seqType='dna')

		self.assertEqual(len(fastqSequences), 5)
		self.assertEqual(fastqSequences[0].seq[:4], 'CCTA')
		self.assertEqual(fastqSequences[0].quality.string[:4], '@8@>')
		self.assertEqual(fastqSequences[0].quality.scores[:4], [31, 23, 31, 29])

class QualityTest(unittest.TestCase):

	def test_quality(self):
		phred33 = seq.QualityScore('!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI', 'phred33')
		phred64 = seq.QualityScore('@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh', 'phred64')

		self.assertEqual(phred33.scores, list(range(41)))
		self.assertEqual(phred64.scores, list(range(41)))

	def test_plot_average_quality_score_per_base(self):
		seq.plot_average_quality_score_per_base('./tests/data/dna.fastq')

if __name__ == '__main__':
	seq.plot_average_quality_score_per_base('./data/dna.fastq')