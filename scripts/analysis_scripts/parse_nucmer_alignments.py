import sys
from collections import defaultdict, Counter
from operator import itemgetter
from itertools import groupby

def in_flagged_region(contig,loci):
	if contig not in flagged_regions: return False
	else:
		if loci in flagged_regions[contig]: return True
		else:
			for i in flagged_regions[contig]:
				if abs(i-loci) <= 3: return True
			else: return False

def indel_ranges(contig,data):
	# print(data)
	for k,g in groupby(enumerate(data), lambda x:x[0]-x[1]):
		x = list(map(itemgetter(1), g))
		L = len(x)
		if L > 20:
			error_types['long_indel'] += 1
			out2.write('%s %d %d %d %s\n' % (contig,x[0],x[-1],L,'short_indel'))
		elif L >= 3:
			out2.write('%s %d %d %d %s\n' % (contig,x[0],x[-1],L,'short_indel'))
			error_types['short_indel'] += 1
		else:
			if L == 1:
				W = results[contig][x[0]]
				LL = len(W)
				if LL > 20:
					out2.write('%s %d %d %d %s\n' % (contig,x[0],x[-1],LL,'long_indel'))
					error_types['long_indel'] += 1
				elif LL >= 3:
					out2.write('%s %d %d %d %s\n' % (contig,x[0],x[-1],LL,'short_indel'))
					error_types['short_indel'] += 1
				else:
					if homopolymer(W[0],contig,x[0]): error_types['homopolymer'] += 1
					# elif gc_rich(W[0]): error_types['gc_rich'] += 1
					elif short_repeat(W[0],contig,x[0]): error_types['short_repeat'] += 1
					else:
						out2.write('%s %d %s %s %s %s\n' % (contig,x[0],W[0][:10],W[0][10],W[0][11:],'unknown'))
						error_types['unknown'] += 1
			else:
				W = results[contig][x[0]]
				if homopolymer(W[0],contig,x[0]): error_types['homopolymer'] += 1
				# elif gc_rich(W[0]): error_types['gc_rich'] += 1
				elif short_repeat(W[0],contig,x[0]): error_types['short_repeat'] += 1
				else:
					out2.write('%s %d %s %s %s %s\n' % (contig,x[0],W[0][:10],W[0][10],W[0][11:],'unknown'))
					error_types['unknown'] += 1

def homopolymer(window2,Contig,Loci):
	if '.' in window2:
		window = window2.replace('.','')
		end = 11
	else:
		window = window2
		end = 12
	for i in range(7,end):
		kmer = window[i:i+3]
		if len(set(kmer)) == 1:
			out2.write('%s %d %s %s %s %s %s\n' % (Contig,Loci,window[:10],window[10],window[11:],kmer,'homopolymer'))
			return True
	return False

def short_repeat(window2,Contig,Loci):
	if '.' in window2:
		window = window2.replace('.','')
		end = 14
	else:
		window = window2
		end = 15
	r = defaultdict(list)
	for i in range(4,end):
		kmer = window[i:i+3]
		r[kmer].append(i)
	for k,v in r.items():
		if len(v) > 1:
			for i in range(len(v)-1):
				for j in range(i+1,len(v)):
					if abs(v[i]-v[j]) == 3:
						out2.write('%s %d %s %s %s %s %s\n' % (Contig,Loci,window2[:10],window2[10],window2[11:],k,'3mer_repeat'))
						return True
	r = defaultdict(list)
	for i in range(6,end-1):
		kmer = window[i:i+2]
		r[kmer].append(i)
	for k,v in r.items():
		if len(v) > 1:
			for i in range(len(v)-1):
				for j in range(i+1,len(v)):
					if abs(v[i]-v[j]) == 2:
						out2.write('%s %d %s %s %s %s %s\n' % (Contig,Loci,window2[:10],window2[10],window2[11:],k,'2mer_repeat'))
						return True
	return False

def gc_rich(window):
	window = window.lower()
	for i in range(4,11):
		kmer = window[i:i+6]
		C = Counter(kmer)
		if (C['g'] + C['c']) >= 4:
			# print(window,'gc')
			return True
	return False

delta = sys.argv[1]
snps = delta.replace('delta','snps')
results_out = delta.replace('delta','results')
error_summary = delta.replace('delta','error_summary')
flagged_regions_file = sys.argv[2]
flagged_regions = defaultdict(list)

for i in open(flagged_regions_file):
	contig,loci,error_type = i.strip().split('\t')
	flagged_regions[contig].append(int(loci))

r = {}

for i in open(delta):
	tmp = i.strip().split(' ')
	if len(tmp) >= 4:
		if i[0] == ">":
			ref,query,refLen,queryLen = tmp[0][1:],tmp[1],int(tmp[2]),int(tmp[3])
			if ref not in r:
				r[ref] = {query:[0 for x in range(refLen)]}
			else:
				if query not in r[ref]:
					r[ref][query] = [0 for x in range(refLen)]
		elif len(tmp) == 7:
			ref_start,ref_end = int(tmp[0]),int(tmp[1])
			for j in range(min(ref_start,ref_end),max(ref_start,ref_end)):
				r[ref][query][j] = 1

aln_mapping = {}
total_aln_len = 0
total_assembly_len = 0

for k,v in r.items():
	for z,w in v.items():
		C = sum(w)
		T = len(w)
		if C/T >= 0.95:
			aln_mapping[k] = {z:[C,T,C/T]}
			total_aln_len += C
			total_assembly_len += T
			# print(k,z,C,T,C/T)


results = {}
error_types = {'unknown':0,'homopolymer':0,'short_repeat':0,'gc_rich':0,'long_indel':0,'short_indel':0}

indels,mismatches = 0,0
plasmid_error = 0
chromosome_error = 0

with open(error_summary,'w') as out2:
	for i in open(snps):
		if 'DIST' in i: continue
		tmp = [x for x in i.strip().split(' ') if x != '']
		if len(tmp) == 14:
			loci_ref,ref,query,loci_query,ref_context,query_context,ref_contig,query_contig = int(tmp[0]),tmp[1],tmp[2],int(tmp[3]),tmp[8],tmp[9],tmp[13].split('\t')[0],tmp[13].split('\t')[1]
			if in_flagged_region(ref_contig,loci_ref): continue # commenting this line out will estimate errors without accounting for flagged regions in pacbio assemblies
			if ref_contig in aln_mapping:
				if query_contig in aln_mapping[ref_contig]:
					if '.' in [ref,query]:
						indels += 1
						if 'ctg.s1.000000F' == ref_contig: chromosome_error += 1
						else: plasmid_error += 1
						if ref_contig not in results:
							results[ref_contig] = defaultdict(list)
						results[ref_contig][loci_ref].append(ref_context)
					else:
						mismatches += 1
						if 'ctg.s1.000000F' == ref_contig: chromosome_error += 1
						else: plasmid_error += 1
						if homopolymer(ref_context,ref_contig,loci_ref): error_types['homopolymer'] += 1
						# elif gc_rich(ref_context): error_types['gc_rich'] += 1
						elif short_repeat(ref_context,ref_contig,loci_ref): error_types['short_repeat'] += 1
						else:
							out2.write('%s %d %s %s %s %s\n' % (ref_contig,loci_ref,ref_context[:10],ref_context[10],ref_context[11:],'unknown'))
							error_types['unknown'] += 1

	for k,v in results.items():
		indel_ranges(k,list(v.keys()))

# print(error_types)

total = indels + mismatches

results = [delta.split('/')[-2].split('.')[0],total,chromosome_error,plasmid_error,mismatches,indels,total_aln_len,total_assembly_len,int(total_assembly_len*total/total_aln_len),int(total_assembly_len*mismatches/total_aln_len),int(total_assembly_len*indels/total_aln_len),error_types['long_indel'],error_types['short_indel'],error_types['homopolymer'],error_types['short_repeat'],error_types['unknown']]
results = '\t'.join([str(x) for x in results])

with open(results_out,'w') as out:
	out.write(results+'\n')
