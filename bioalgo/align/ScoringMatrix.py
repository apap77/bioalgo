class ScoringMatrix:
	def __init__(self, matrixId):

		if matrixId in ['blosum62', 'pam250']:
			self.characters = 'ACDEFGHIKLMNPQRSTVWY'
		else:
			self.characters = 'ATCG'

		self.scoreMat = self._parse_scoring_matrix(matrixId)
	
	def __getitem__(self, i):
		return self.scoreMat[i]

	def _parse_scoring_matrix(self, matrixId):
		scoreMat = {}

		with open('bioalgo/align/data/%s.txt' % matrixId) as inFile:
			for c, line in zip(self.characters, inFile.readlines()):
				scoreMat[c] = dict(zip(self.characters, list(map(int, line.strip().split()))))

		return scoreMat

if __name__ == '__main__':
	print(ScoringMatrix('blosum62')['A']['A'])