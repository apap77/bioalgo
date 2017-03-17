import matplotlib.pyplot as plt
from collections import defaultdict

def kmer_dict(seq, k):
	kmerDict = defaultdict(list)
	for i in range(len(seq) - k + 1):
		kmerDict[seq[i:i+k]].append(i)

	return kmerDict

def reverse_complement(seq):
	complement = str.maketrans('ACGT', 'TGCA')
	return seq.translate(complement)[::-1]

def dot_plot(seq1, seq2, k=10):
	kmerDict = kmer_dict(seq1, k)

	dots = []
	for i in range(len(seq2) - k + 1):
		if seq2[i:i+k] in kmerDict:
			for index in kmerDict[seq2[i:i+k]]:
				dots.append((index, i))
		if reverse_complement(seq2[i:i+k]) in kmerDict:
			for index in kmerDict[reverse_complement(seq2[i:i+k])]:
				dots.append((index, i))

	x, y = zip(*dots)

	plt.scatter(x, y)
	plt.show()


if __name__ == '__main__':
	seq1 = 'AGCAGGTTATCTTCCTGT'
	seq2 = 'AGCAGGTTATCTTCCTGT'

	dot_plot(seq1, seq2, k=3)
