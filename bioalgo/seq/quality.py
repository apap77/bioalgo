class QualityScore:

	def __init__(self, string, encoding='phred33'):
		self.string = string
		self.encoding = encoding
		self._parse_string()

	def _parse_string(self):
		if self.encoding == 'phred33':
			self.scores = [ord(c) - 33 for c in self.string]

		elif self.encoding == 'phred64':
			self.scores = [ord(c) - 64 for c in self.string]

def plot_average_quality_score_per_base(filePath):
	import matplotlib.pyplot as plt
	import numpy as np
	from bioalgo.seq.parse import Fastq

	sequences = Fastq(filePath)

	# compute average quality score per base
	scores = np.array([sequence.quality.scores for sequence in sequences])
	meanScores = np.mean(scores, axis=0)

	# show plot
	plt.plot(meanScores)
	plt.title('Average quality score per base')
	plt.xlabel('Base Position')
	plt.ylabel('Phred Score')
	plt.show()
	