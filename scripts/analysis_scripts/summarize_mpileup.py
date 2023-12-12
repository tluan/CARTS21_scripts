import sys
from collections import Counter
from scipy.stats import entropy
import numpy as np
import re

def parse_calls(calls):
	calls = calls.replace('$','').replace('^','')
	insertions,deletions = calls.count('+'),calls.count('-')
	calls = re.split('\+|-', calls)
	new = ''
	for i in calls:
		cc = '0'
		for j in i:
			try:
				int(cc+j)
				cc += j
			except ValueError:
				break
		cc = int(cc)
		new += i[cc+1:]
	return new,insertions,deletions

'''
def parse_calls(calls):
	nums = sorted([int(x) for x in re.findall(r'\d+', calls)],reverse=True)
	if nums == []: return calls
	for x in nums:
		x = str(x)
		calls = calls[:calls.find(x)]+calls[calls.find(x)+int(x)+1:]
	return calls
'''

mpileup = sys.argv[1]

with open(mpileup+'.variable_loci.txt','w') as out:
	out.write("Contig\tLoci\tTotal_reads\tReads_A\tReads_C\tReads_G\tReads_T\tReads_ins\tReads_del\tNum_variants\tEntropy\tbaseQual\tmapQual\n")
	for h,i in enumerate(open(mpileup)):
		tmp = i.strip().split('\t')
		contig,loci,doc,calls,baseQuals,mapQuals = tmp[0],tmp[1],int(tmp[3]),tmp[4].lower(),tmp[5],tmp[6]
		baseQuals = np.mean([ord(x)-33 for x in baseQuals])
		mapQuals = np.mean([ord(x)-33 for x in mapQuals])
		calls,insertions,deletions = parse_calls(calls)
		C = Counter(calls)
		r = []
		V = 0
		total = 0
		for j in 'acgt':
			c = C.get(j,0)
			total += c
			r.append(str(c))
			if c > 0: V += 1
		r += [str(insertions), str(deletions)]
		total += (insertions + deletions)
		if total == 0: shannon_entropy = 0
		else:
			r2 = [float(x)/total for x in r]
			shannon_entropy = entropy(r2,base=2)
		out.write("%s\t%s\t%d\t%s\t%d\t%f\t%d\t%d\n" % (contig,loci,total,'\t'.join(r),V,shannon_entropy,baseQuals,mapQuals))
