class ScoringMatrix:
	def __init__(self, matrixId):
		self.aminoAcids = "ACDEFGHIKLMNPQRSTVWY"
		self.scoreMat = self._parse_scoring_matrix(matrixId)
	
	def __getitem__(self, i):
		return self.scoreMat[i]

	def _parse_scoring_matrix(self, matrixId):
		scoreMat = {}
		with open('./data/%s.txt' % matrixId) as inFile:
			for aa, line in zip(self.aminoAcids, inFile.readlines()):
				scoreMat[aa] = dict(zip(self.aminoAcids, list(map(int, line.strip().split()))))

		return scoreMat

if __name__ == '__main__':
	print(ScoringMatrix('blosum62')['A']['A'])