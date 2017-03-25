from collections import defaultdict

def all_kmers_of(text, k):
	for i in range(len(text) - k + 1):
		yield text[i:i+k]

class DeBrujin:

	def __init__(self, text, k):
		self.graph = defaultdict(list)

		for kmer in all_kmers_of(text, k):
			self.graph[kmer[:-1]].append(kmer[1:])

	def __getitem__(self, node):
		return self.graph[node]

	def __str__(self):
		return '\n'.join(['%s: %s' % (k, ', '.join(v)) for k, v in sorted(self.graph.items())])