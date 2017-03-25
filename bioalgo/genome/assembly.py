from collections import defaultdict
from bioalgo.genome.util import all_kmers_of, eulerian_path

class DeBrujin:

	def __init__(self, text, k=3):
		self.graph = defaultdict(list)

		if type(text) is str:
			self._construct_graph_from_text(text, k)
		if type(text) is list :
			self._construct_graph_from_patterns(text, k)

		self.edges = []
		for key, values in self.graph.items():
			self.edges.extend([(key, value) for value in values])

	def reconstruct(self):
		eulerianPath = eulerian_path(self.edges)
		return eulerianPath[0][0] + ''.join([src[-1] for src, dst in eulerianPath[1:]]) + eulerianPath[-1][1][-1]

	def _construct_graph_from_text(self, text, k):
		for kmer in all_kmers_of(text, k):
			self.graph[kmer[:-1]].append(kmer[1:])

	def _construct_graph_from_patterns(self, patterns, k):
		for pattern in patterns:
			for kmer in all_kmers_of(pattern, k):
				self.graph[kmer[:-1]].append(kmer[1:])

	def __getitem__(self, node):
		return self.graph[node]

	def __str__(self):
		return '\n'.join(['%s -> %s' % (key, ', '.join(values)) for key, values in sorted(self.graph.items())])