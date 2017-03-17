import sys; sys.path.append('../')
import bioalgo.seq as seq

if __name__ == '__main__':
	d = seq.RNA('AUGAAUGAUUGA')
	kc = d.kmer_composition(k=3)
	print(kc)
	print(d.translate(frameStart=-1))
	print(d.reverse_complement())

	d = seq.DNA('ATGTTATGGATAG')
	kc = d.kmer_composition(k=2)
	print(kc)