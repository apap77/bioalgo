import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm

def kmer_dict(seq, k):
	kmerDict = defaultdict(list)
	for i in tqdm(range(len(seq) - k + 1)):
		kmerDict[seq[i:i+k]].append(i)

	return kmerDict

def reverse_complement(seq):
	complement = str.maketrans('ACGTN', 'TGCAN')
	return seq.translate(complement)[::-1]

def dot_plot(seq1, seq2, k=10):
	kmerDict = kmer_dict(seq1, k)

	dots = set()
	for i in tqdm(range(len(seq2) - k + 1)):
		for index in kmerDict[seq2[i:i+k]] + kmerDict[reverse_complement(seq2[i:i+k])]:
			dots.add((index, i))

	x, y = zip(*list(dots))

	plt.scatter(x, y, marker='.', color='black')
	plt.show()


if __name__ == '__main__':
	seq1 = 'AGCAGGTTATCTTCCTGT'
	seq2 = 'AGCAGGTTATCTTCCTGT'

	dot_plot(seq1, seq2, k=3)
