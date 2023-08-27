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
parser.add_argument("-i", "--i", help = "String; Name of input file containing junctions in bed format.")
parser.add_argument("-o", "--o", help = "String; Path to output directory to place output quant files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-b", "--b", help = "String; Path to directory with input bam files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-bams", "--bams", help = "List: List of bam filename prefixs.", nargs = "+")
parser.add_argument("-s", "--s", help = "String: suffix (excluding .bam) for bam filenames. Default is ''", default = "")
parser.add_argument("-l", "--l", help = "String: Library orientation. Default is u.", default = "u")

args = parser.parse_args()
"""
chr	strand	ret.intron.start	ret.intron.end	unqiue.regions	ovlap.exons
0	   1		2		3			4		5
1	1148473	1149042	D2.event739.1	0		-	
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

def check_read_compatible(bamObj, r, don, acc, strand, libtype):
	"""
	read = Aligned segment object of pysam that needs to be checked for correct orientation
	strand = string; can be + or - depending on strand of region of interest in which reads are being counted
	libtype = string representing strandedness of your RNA seq library
	intron = tuple of intron coordinates (junction coordinates) to count the split reads mapped across the intron (junction)
	Value: Returns True or False depending on whether the read fullfills all criteria
	"""
	if (strand_consistent(r, strand, libtype)) and (r.cigarstring.count("N") >= 1):
		for irn in bamObj.find_introns([r]):
			if (don == irn[0]) and (irn[1] == acc):
				return True

	return False

def count_mmAdj(bams, ridx, libtype):
	bam_wise_counts = []

	for bp in bams:
		with pys.AlignmentFile(bp, "rb") as bamObj:
			chrom, strand, don, acc = ridx.split("__") 
					
			reads = [r for r in bamObj.fetch(chrom, int(don) - 5, int(acc) + 5) if check_read_compatible(bamObj, r, int(don), int(acc), strand, libtype)]
			num_reads = len(reads) 
				
			rnames = [r.query_name for r in reads]
			counts = Counter(rnames)
						
			nh = [r.get_tag("NH") for r in reads]
			ncountAdj = np.round(sum([1.0/(counts[rnames[k]]*nh[k]) for k in range(num_reads)]), 2)

			bam_wise_counts.append(ncountAdj)

	return bam_wise_counts

jpath = os.path.join(args.w, args.i)
bams = [os.path.join(args.b, bam+args.s+".bam") for bam in args.bams]
libtype = args.l

bamcols = []
for b in args.bams:
	bamcols.extend([b+".count"])

opf = os.path.join(args.o, args.i+".countsAdj.out")

with open(opf, "w") as oph:
	oph.write("\t".join(["chr","jn.start", "jn.end", "strand", "id"] + bamcols) + "\n")
	with open(jpath, "r") as reth:
		for line in reth:
			rec = line.strip().split("\t")

			if rec[0] == "chr":
				continue
			
			jnidx = "__".join([rec[0], rec[5], rec[1], rec[2]])
			chrom, jstart, jend, strand = rec[0], rec[1], rec[2], rec[5]

			jncounts = count_mmAdj(bams, jnidx, libtype)
			
			orec = [chrom, strand, jstart, jend, jnidx] + [str(c) for c in jncounts]
			oph.write("\t".join(orec)+"\n")
			
