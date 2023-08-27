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
parser.add_argument("-l", "--l", help = "String: Library orientation. Default is u (un-stranded).", default = "u")

args = parser.parse_args()


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
			if (don == irn[0]) or (irn[1] == acc):
				return True

	return False

def count_mmAdj(bams, ridx, libtype):
	bam_wise_counts = []

	for bp in bams:
		with pys.AlignmentFile(bp, "rb") as bamObj:
			chrom, strand, don, acc = ridx.split("__")
			reads = None
			try:
				reads = [r for r in bamObj.fetch(chrom, int(don) - 5, int(acc) + 5) if check_read_compatible(bamObj, r, int(don), int(acc), strand, libtype)]
			except:
				reads = [r for r in bamObj.fetch("chr"+chrom, int(don) - 5, int(acc) + 5) if check_read_compatible(bamObj, r, int(don), int(acc), strand, libtype)]
			num_reads = len(reads) 
				
			rnames = [r.query_name for r in reads]
			counts = Counter(rnames)
						
			nh = [r.get_tag("NH") for r in reads]
			ncountAdj = np.round(sum([1.0/(counts[rnames[k]]*nh[k]) for k in range(num_reads)]), 2)

			bam_wise_counts.append(ncountAdj)

	return bam_wise_counts


retjpath = os.path.join(args.w, args.i)
retJns = {} # initialize empty retjns dictionary >> needs to be filled as -- retJns[idx] = {"tx":[], "ribases":[], "unq":[]}, idx = "__".join([chr, strand, start, end]) ## chr, strand, start (0-b), end (1-b)

with open(retjpath, "r") as reth:
	for line in reth:
		rec = line.strip().split("\t")

		if rec[0] == "chr":
			continue
		
		idx = "__".join(rec[0:4])
		tx = rec[4]
		unq_regions = rec[5]

		if idx not in retJns:
			retJns[idx] = {"tx":[], "unq":[]}

		retJns[idx]["tx"].append(tx)
		retJns[idx]["unq"] = unq_regions

## get the skipping counts across the junctions
skipcols = []
for b in args.bams:
	skipcols.extend([b+".skip_counts"])

libtype = args.l

bams = [os.path.join(args.b, bam+args.s+".bam") for bam in args.bams]
with open(os.path.join(args.o, "SkippingCounts.out"), "w") as skiph:
	
	skiph.write("\t".join(["chr", "strand", "ret.intron.start", "ret.intron.end"] + skipcols) + "\n")
	
	for ridx in retJns:
		chrom, strand, istart, iend = ridx.split("__")
		skip_counts = count_mmAdj(bams, ridx, libtype)
		orec = [chrom, strand, istart, iend] + [str(c) for c in skip_counts]
		skiph.write("\t".join(orec)+"\n")

