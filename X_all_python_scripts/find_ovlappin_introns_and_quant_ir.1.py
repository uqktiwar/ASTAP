import pysam as pys
import sys
import numpy as np
import statistics 
from collections import Counter 
import argparse as ap
import os
from operator import itemgetter

parser = ap.ArgumentParser()
parser.add_argument("-w", "--w", help = "String; Path to working directory with input file. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-i", "--i", help = "String; Name of input file contaiining retained junctions with their unique regions.")
parser.add_argument("-o", "--o", help = "String; Path to output directory to place output quant files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-b", "--b", help = "String; Path to directory with input bam files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-bams", "--bams", help = "List: List of bam filename prefixs.", nargs = "+")
parser.add_argument("-s", "--s", help = "String: suffix (excluding .bam) for bam filenames. Default is ''", default = "")
parser.add_argument("-l", "--l", help = "String: Library orientation. Default is rf.", default = "rf")

args = parser.parse_args()
"""
chr	strand	ret.intron.start	ret.intron.end	txids	unqiue.regions	ovlap.exons
0	1		2		3		4				5					6
1	-	1148473	1149042	D2.event739.1	[1148473, 1149042]_[154143187, 154143888]_[154143964, 154144357]	
1	-	1670692	1670929	D2.ovevent113.1	[1670692, 1670929]	
1	-	1683836	1683880	D2.event760.1	[1683836, 1683880]	
"""

def strand_consistent(read, strand, libtype):
	# checking reads according to the rna seq library strandedness
	if libtype == "rf":
		if strand == "+":
			if read.is_read1 and not read.is_reverse:
				return False
			if read.is_read2 and read.is_reverse:
				return False
		if strand == "-":
			if read.is_read1 and read.is_reverse:
				return False
			if read.is_read2 and not read.is_reverse:
				return False
		
	if libtype == "fr":
		if strand == "+":
			if read.is_read2 and not read.is_reverse:
				return False
			if read.is_read1 and read.is_reverse:
				return False
		if strand == "-":
			if read.is_read2 and read.is_reverse:
				return False
			if read.is_read1 and not read.is_reverse:
				return False
	return True


retjpath = os.path.join(args.w, args.i)
retJns = {} # initialize empty retjns dictionary >> needs to be filled as -- retJns[idx] = {"tx":[], "ribases":[], "unq":[]}, idx = "__".join([chr, strand, start, end]) ## chr, strand, start (0-b), end (1-b)

with open(retjpath, "r") as reth:
	for line in reth:
		rec = line.strip().split("\t")

		if rec[0] == "chr":
			continue
		
		idx = "__".join(rec[0:4])
		txs = rec[4]
		
		unq_regions = rec[5]
		
		for tx in txs.split(","):
			if tx not in retJns:
				retJns[tx] = {}
			
			if idx not in retJns[tx]:
				retJns[tx][idx] = []

			retJns[tx][idx] = unq_regions


## group introns by their unq regions
group_by_unq = {}
all_unq_idxs = {}
for tx in retJns:
	if tx not in group_by_unq:
		group_by_unq[tx] = {}
	
	if tx not in all_unq_idxs:
		all_unq_idxs[tx] = []

	for ri in retJns[tx]:
		chrom, strand, ris, rie = ri.split("__")

		if retJns[tx][ri] == "No_Unq_Regions":
			continue
		
		unq_idx = "__".join([chrom, strand, retJns[tx][ri]])
		unqs = [[int(ue) for ue in uq.replace("[","").replace("]","").split(", ")] for uq in retJns[tx][ri].split("_")]
		span_unq = max([uq[1] for uq in unqs]) - min([uq[0] for uq in unqs])
		ri_unq = int(rie) - int(ris)

		all_unq_idxs[tx].append([unq_idx, ri, span_unq, ri_unq])
		
		if unq_idx not in group_by_unq[tx]:
			group_by_unq[tx][unq_idx] = {"counts":[], "quantiles":[], "stdev":[]}


## remove introns whose unq_idx are a subset of the unq_idxs of other introns within the same transcript
final_groups = {}
for tx in all_unq_idxs:
	final_groups[tx] = {}

	all_unq_idxs[tx] = sorted(all_unq_idxs[tx], key = itemgetter(2, 3), reverse = True)
	rmkeys = []
	for i in range(len(all_unq_idxs[tx])-1):
		ut1 = all_unq_idxs[tx][i]
		#print(ut1)
		uidx1 = ut1[0]
		ri1 = ut1[1]
		
		for j in range(i+1, len(all_unq_idxs[tx])):
			ut2 = all_unq_idxs[tx][j]
			uidx2 = ut2[0]
			ri2 = ut2[1]

			if (uidx1 == uidx2) and (ut1[1] == ut2[1]):
				continue
					
			setu1, setu2 = set(uidx1.split("__")[-1].split("_")), set(uidx2.split("__")[-1].split("_"))
						
			if setu2.issubset(setu1):
				if setu1.issubset(setu2):
					retJns[tx].pop(ri2, None)
					continue
				rmkeys.append(uidx2)

	for uidx in group_by_unq[tx]:
		if uidx not in rmkeys:
			final_groups[tx][uidx] = group_by_unq[tx][uidx]


final_final_groups = {}

for tx in final_groups:
	for unq_idx in final_groups[tx]:
		if unq_idx not in final_final_groups:
			final_final_groups[unq_idx] = {"counts":[], "quantiles":[], "stdev":[]}




bams = [os.path.join(args.b, bam+args.s+".bam") for bam in args.bams]
for b in bams:
	samfile = pys.AlignmentFile(b, "rb")
	for uidx in final_final_groups:
		chrom, strand, unqlist = uidx.split("__")
			
		final_counts_list = []
		for tunq in unqlist.split("_"):
				##print(tunq)
			tunq = [int(e) for e in tunq.replace("[", "").replace("]", "").replace(" ", "").split(",")]
			a, b = tunq
			b = max(b, a+1)
			pos_wise_adj_counts = [0.0]*(b - a)
			pos_to_update = list(range(a, b))

			pcoliter = None
			try:
				pcoliter = samfile.pileup(chrom, a, b, truncate = True)
			except:
				pcoliter = samfile.pileup("chr"+chrom, a, b, truncate = True)
				
				#print(iter)
			for pcol in pcoliter:
				reads = []
				ftcount = 0.0
				for pread in pcol.pileups:
					if (not pread.is_del and not pread.is_refskip) and (strand_consistent(pread.alignment, strand, libtype = args.l)):
						reads.append(pread.alignment)
					
				if reads != []:
					num_reads = len(reads)
					rnames = [r.query_name for r in reads]
					counts = Counter(rnames)
								
					nh = [r.get_tag("NH") for r in reads]
					ftcount = np.round(sum([1.0/(counts[rnames[k]]*nh[k]) for k in range(num_reads)]), 2)

				pos_idx = pos_to_update.index(pcol.reference_pos) ## 0 based
				pos_wise_adj_counts[pos_idx] = ftcount
				
				##print(len(pos_wise_adj_counts), b-a)
			final_counts_list.extend(pos_wise_adj_counts)
				
		med_count = statistics.median(final_counts_list)
		quantiles = (np.quantile(final_counts_list, 0.25), np.quantile(final_counts_list, 0.5), np.quantile(final_counts_list, 0.75))
			##print(np.sum(final_counts_list), len(final_counts_list))
		std = np.std(final_counts_list)
			
		final_final_groups[uidx]["counts"].append(med_count) 
		final_final_groups[uidx]["quantiles"].append(quantiles)
		final_final_groups[uidx]["stdev"].append(std)
		
	samfile.close()


bamcols = []
for b in args.bams:
	bamcols.extend([b+".count", b+".quantl", b+".stdev"])


opf = os.path.join(args.o, args.i+".ri.countsAdj.out")


with open(opf, "w") as oph:
	oph.write("\t".join(["chr", "strand", "ret.intron.start", "ret.intron.end", "isoforms", "unqiue.regions"] + bamcols) + "\n")
	for tx in retJns:
		for ri in retJns[tx]:
			chrom, strand, ris, rie = ri.split("__")
			unq_regions = retJns[tx][ri]
			
			if unq_regions == "No_Unq_Regions":
				oph.write("\t".join([chrom, strand, ris, rie, tx, unq_regions] + ["NA"]*len(bamcols)) + "\n")
				continue

			uidx = "__".join([chrom, strand, unq_regions])
			if uidx in final_final_groups:
				counts = final_final_groups[uidx]["counts"]
				quants = final_final_groups[uidx]["quantiles"]
				stdev = final_final_groups[uidx]["stdev"]

				ocounts = []
				for i in range(len(args.bams)):
					ocounts.extend([str(e) for e in [counts[i], quants[i], stdev[i]]])

				oph.write("\t".join([chrom, strand, ris, rie, tx, unq_regions] + ocounts) + "\n")




		



