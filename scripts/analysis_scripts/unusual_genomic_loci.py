import sys
import numpy as np
import itertools
from collections import defaultdict

def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

stats = sys.argv[1] # variable loci report from mpileup analysis
min_window_size = 5

results = {}
r_mq = defaultdict(list)
r_doc = defaultdict(list)

with open(sys.argv[1]+'.flagged_regions','w') as out:
	for h,i in enumerate(open(stats)):
		if h ==	0: continue
		tmp = i.strip().split('\t')
		contig,loci,doc,A,C,G,T,insertion,deletion,ent,mapq = tmp[0],int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]),int(tmp[6]),int(tmp[7]),int(tmp[8]),float(tmp[10]),int(tmp[-1])
		if contig not in results: results[contig] = {}
		indel = insertion+deletion
		if mapq < 40:
			r_mq[contig].append(loci)
		if doc < 20:
			r_doc[contig].append(loci)
		M = max([A,C,G,T,insertion,deletion])
		if M > 0:
			c = 0
			for x in [A,C,G,T,insertion,deletion]:
				if x/M > 0.4:
					c += 1
			if c > 1:
				if loci not in results[contig]: results[contig][loci] = ["potential SNP"]
				else: results[contig][loci].append("potential SNP")
	
	for k,v in r_mq.items():
		if k not in results: results[k] = {}
		for i in intervals_extract(v):
			d = abs(i[0]-i[1])
			if d >= min_window_size:
				for j in range(min(i[0],i[1]),max(i[0],i[1])):
					if j not in results[k]: results[k][j] = ["low mapq score"]
					else: results[k][j].append("low mapq score")
	
	for k,v in r_doc.items():
		if k not in results: results[k] = {}
		for i in intervals_extract(v):
			d = abs(i[0]-i[1])
			if d >= min_window_size:
				for j in range(min(i[0],i[1]),max(i[0],i[1])):
					if j not in results[k]: results[k][j] = ["low DOC"]
					else: results[k][j].append("low DOC")
	
	for contig,v in sorted(results.items()):
		for z,w in v.items():
			out.write("%s\t%d\t%s\n" % (contig,z,",".join(w)))


