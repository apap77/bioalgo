from bioalgo.align.ScoringMatrix import ScoringMatrix

class BasicGlobalAligner:

	def __init__(self, seq1, seq2, scoringMatrix='blosum62', indelPenalty=-1):
		self.seq1 = seq1
		self.seq2 = seq2
		self.indelPenalty = indelPenalty
		self.scoringMatrix = ScoringMatrix(scoringMatrix)


	def align(self):
		self._initialize_alignment_matrix()
		self._initialize_backtrack_matrix()

		self._run_alignment_algorithm()
		self._backtrack()


	def _initialize_alignment_matrix(self):
		self.alignmentMatrix = [[0] * (len(self.seq2) + 1) for _ in range(len(self.seq1) + 1)]

		for i in range(1, len(self.seq1) + 1):
			self.alignmentMatrix[i][0] = self.indelPenalty * i

		for j in range(1, len(self.seq2) + 1):
			self.alignmentMatrix[0][j] = self.indelPenalty * j


	def _initialize_backtrack_matrix(self):
		self.backtrackMatrix = [[-1] * (len(self.seq2) + 1) for _ in range(len(self.seq1) + 1)]

		# 0 : UP, 1 : LEFT, 2 : UPPER-LEFT in backtrack matrix
		for i in range(1, len(self.seq1) + 1):
			self.backtrackMatrix[i][0] = 0

		for j in range(1, len(self.seq2) + 1):
			self.backtrackMatrix[0][j] = 1


	def _run_alignment_algorithm(self):
		for i in range(1, len(self.seq1)+1):
			for j in range(1, len(self.seq2)+1):
				backtrackValue, alignmentScore = self._determine_cell_entry(i, j)

				self.alignmentMatrix[i][j] = alignmentScore
				self.backtrackMatrix[i][j] = backtrackValue

		self.score = self.alignmentMatrix[len(self.seq1)][len(self.seq2)]


	def _determine_cell_entry(self, i, j):
		scoreFromUpperCell = self.alignmentMatrix[i-1][j] + self.indelPenalty
		scoreFromLeftCell = self.alignmentMatrix[i][j-1] + self.indelPenalty
		scoreFromUpperLeftCell = self.alignmentMatrix[i-1][j-1] + self.scoringMatrix[self.seq1[i-1]][self.seq2[j-1]]

		candidates = [scoreFromUpperCell, scoreFromLeftCell, scoreFromUpperLeftCell]
		return max([(backtrackValue, alignmentScore) for backtrackValue, alignmentScore in enumerate(candidates)], key=lambda x: x[1])


	def _backtrack(self):
		i, j = len(self.seq1), len(self.seq2)
		seq1Augmentation, seq2Augmentation = [], []

		while i != 0 or j != 0:
			if self.backtrackMatrix[i][j] == 0:
				seq1Augmentation.append(self.seq1[i-1])
				seq2Augmentation.append('-')
				i, j = i-1, j
			elif self.backtrackMatrix[i][j] == 1:
				seq1Augmentation.append('-')
				seq2Augmentation.append(self.seq2[j-1])
				i, j = i, j-1
			else:
				seq1Augmentation.append(self.seq1[i-1])
				seq2Augmentation.append(self.seq2[j-1])
				i, j = i-1, j-1

		self.seq1Result = ''.join(seq1Augmentation[::-1])
		self.seq2Result = ''.join(seq2Augmentation[::-1])

	def __str__(self):
		try:
			return '\n'.join([self.seq1Result, self.seq2Result])
		except AttributeError:
			raise ValueError('Please call align() first.')