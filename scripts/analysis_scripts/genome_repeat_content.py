import numpy as np
from Bio import SeqIO

def find_repeats(seq,klen):
	if klen in [2,3]:
		min_len = klen*2
		min_set_len = 2
	elif klen == 1:
		min_len = 3
		min_set_len = 1
	repeats = []
	i = 0
	while i < len(seq) - klen:
		if i in exclude:
			i += 1
			continue
		tmp = seq[i:i+klen]
		for j in range(klen,100000,klen):
			if seq[i:i+klen] == seq[i+j:i+j+klen]:
				tmp += seq[i:i+klen]
			else:
				if len(tmp) >= min_len:
					if len(set(tmp)) >= min_set_len:
						repeats.append([tmp,j])
						exclude.update([i])
					i += j+klen
				else:
					i += 1
				break
	return repeats

fastas = 'pacbio_assembly_files.txt'
samples = [fasta.strip() for fasta in open(fastas)]

for fasta in samples:
	exclude = set()
	for klen in range(1,4):
		tmp_results = []
		for i in SeqIO.parse(fasta,'fasta'):
			s = str(i.seq)
			tmp_results += find_repeats(s,klen)
		distribution_repeats = [len(x[0]) for x in tmp_results]
		print(fasta.split('/')[-1],klen,len(tmp_results),sum(distribution_repeats),min(distribution_repeats),np.median(distribution_repeats),max(distribution_repeats))
