from bioalgo.seq.codon import CodonTable
from collections import Counter

class Seq:
	def __init__(self, seq, description, quality=None):
		self.seq = seq
		self.description = description
		self.quality = quality
		self.elements = None

	def _validate_element(self):
		n, s = set(self.elements), set(self.seq)
		assert s.issubset(n), 'Invalid element %s detected. Please check.' % ', '.join(s - n)

	def kmer_composition(self, k):
		counter = Counter()
		for i in range(len(self.seq) - k + 1):
			counter[self.seq[i:i+k]] += 1

		return counter

	def __str__(self):
		return self.seq

class DNA(Seq):
	def __init__(self, seq, description=None, quality=None):
		super().__init__(seq, description, quality)
		self.elements = ['A', 'T', 'G', 'C']

		self._validate_element()

		self.dnaComplement = str.maketrans('ATGC', 'TACG')
		self.rnaComplement = str.maketrans('ATGC', 'UACG')

	def reverse_complement(self):
		return self.seq.translate(self.dnaComplement)[::-1]

	def transcribe(self, reverse=True):
		complement = self.seq.translate(self.rnaComplement)

		return complement[::-1] if reverse else complement
		
class RNA(Seq):
	def __init__(self, seq, description=None, quality=None):
		super().__init__(seq, description, quality)
		self.elements = ['A', 'U', 'G', 'C']

		self._validate_element()

		self.dnaComplement = str.maketrans('AUGC', 'TACG')
		self.rnaComplement = str.maketrans('AUGC', 'UACG')

	def reverse_complement(self):
		return self.seq.translate(self.rnaComplement)[::-1]

	def translate(self, frameStart=1):
		table = CodonTable()
		s = self.seq if frameStart > 0 else self.seq[::-1]
		aas = [table[s[i:i+3]] for i in range(abs(frameStart-1), len(s), 3) if i+3 <= len(s)]

		return ''.join(aas)

class Peptide(Seq):
	def __init__(self, seq, description=None, quality=None):
		super().__init__(seq, description, quality)
		self.elements = ['A', 'U', 'G', 'C']

		self._validate_element()

if __name__ == '__main__':
	d = DNA('ATTG')

	print(d.transcribe())

