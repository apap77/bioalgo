from bioalgo.seq.quality import QualityScore
from bioalgo.seq.seq import DNA, RNA, Peptide

def Fasta(filePath, seqType='dna'):
	Seq = determine_sequence_type(seqType)

	with open(filePath) as inFile:
		header = inFile.readline().strip()
		assert header.startswith('>'), 'Invalid fasta file.'

		seqs, lines = [], []
		line = inFile.readline().strip()
		while True:
			if line.startswith('>'):
				seqs.append(Seq(''.join(lines), description=header))
				header = line
				lines = []

			elif not line:
				seqs.append(Seq(''.join(lines), description=header))
				header = line
				lines = []
				break
				
			else:
				lines.append(line)
			line = inFile.readline().strip()

	return seqs

def Fastq(filePath, seqType='dna', qualityEncoding='phred33'):
	Seq = determine_sequence_type(seqType)

	seqs = []
	with open(filePath) as inFile:
		while True:
			seq = read_single_fastq_sequence(inFile, Seq, qualityEncoding)
			if not seq:
				break

			seqs.append(seq)		

	return seqs

def determine_sequence_type(seqType):
	if seqType == 'rna':
		Seq = RNA
	elif seqType == 'peptide':
		Seq = Peptide
	else:
		Seq = DNA

	return Seq

def read_single_fastq_sequence(inFile, Seq, qualityEncoding):
	header = inFile.readline().strip()
	if header == '':
		return None

	assert header.startswith('@'), 'Invalid fastq file.'
	seqString = inFile.readline().strip()

	assert inFile.readline().startswith('+'), 'Invalid fastq file.'
	qualityString = inFile.readline().strip()
	quality = QualityScore(qualityString, qualityEncoding)

	return Seq(seq=seqString, description=header, quality=quality)