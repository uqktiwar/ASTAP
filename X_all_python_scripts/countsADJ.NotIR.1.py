import pysam as pys
import argparse as ap
import os
import numpy as np
from collections import Counter

parser = ap.ArgumentParser()
parser.add_argument("-w", "--w", help = "String; Path to working directory with input bed file. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-i", "--i", help = "String; Name of input bed file.")
parser.add_argument("-o", "--o", help = "String; Path to output directory to place output bam files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-b", "--b", help = "String; Path to directory with input bam files. Default is os.getcwd()", default = os.getcwd())
parser.add_argument("-bams", "--bams", help = "List: List of bam filename prefixs.", nargs = "+")
parser.add_argument("-s", "--s", help = "String: suffix (excluding .bam) for bam filenames. Default is ''", default = "")
parser.add_argument("-d", "--d", help = "Int: Degree of events. Default is 2.", default = "2")
parser.add_argument("-l", "--l", help = "String: Library orientation. Default is rf.", default = "rf")

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


def proc_bed(bedPath):
	events = {}

	with open(bedPath) as bedh:
		for line in bedh:
			rec = line.strip().split("\t")
			eid = ".".join(rec[3].split(".")[:-1])
			evs, eve = int(rec[1]), int(rec[2])

			if eid not in events:
				events[eid] = {"chrom":rec[0], "strand":rec[5], "minD":np.inf, "maxA":0 ,"iso":{}}
			
			exons = [(evs+int(s), evs+int(s)+int(l)) for (s,l) in zip(rec[11].split(",")[:-1], rec[10].split(",")[:-1])]
			introns = [(exons[i-1][1], exons[i][0]) for i in range(1, len(exons))]

			try:
				events[eid]["minD"] = min(events[eid]["minD"], introns[0][0])
				events[eid]["maxA"] = max(events[eid]["maxA"], introns[-1][1])
			except:
				pass
			
			events[eid]["iso"][rec[3]] = [introns, []] 
	
	return events

def check_read_compatible(bamObj, r, op, don, acc, introns, strand, libtype):
	"""
	read = Aligned segment object of pysam that needs to be checked for correct orientation
	strand = string; can be + or - depending on strand of region of interest in which reads are being counted
	libtype = string representing strandedness of your RNA seq library
	intron = tuple of intron coordinates (junction coordinates) to count the split reads mapped across the intron (junction)
	Value: Returns True or False depending on whether the read fullfills all criteria
	"""
	if (not strand_consistent(r, strand, libtype)) or (r.cigarstring.count("N") < 1):
		return op
		
	for irn in bamObj.find_introns([r]):
		if (don <= irn[0] < acc) or (don < irn[1] <= acc):
			if irn not in introns:
				return op
	
	for irn in bamObj.find_introns([r]):
		if irn in introns:
			###idx = introns.index(irn)
			op[irn].append(r)

	return op


bedpath = os.path.join(args.w, args.i)
opf = os.path.join(args.o, args.i+".countsAdj.out")
print(opf)
events = proc_bed(bedpath)

bams = [os.path.join(args.b, bam+args.s+".bam") for bam in args.bams]

##print(bams)
##print(obams)
totalEvents = len(events.keys())
eventsFcounts = {}
for bp in bams:
	with pys.AlignmentFile(bp, "rb") as bamObj:
		for eid in events:
			#print(eid)
			chrom, strand = events[eid]["chrom"], events[eid]["strand"]
			don, acc = events[eid]["minD"], events[eid]["maxA"]
			
			if eid not in eventsFcounts:
				eventsFcounts[eid] = {"isos":{}}

			for iso in events[eid]["iso"]:
				introns = events[eid]["iso"][iso][0]
				
				if iso not in eventsFcounts[eid]["isos"]:
					eventsFcounts[eid]["isos"][iso] = {}

				op = {}
				for irn in introns:
					if irn not in eventsFcounts[eid]["isos"][iso]:
						eventsFcounts[eid]["isos"][iso][irn] = []
					op[irn] = []
				
				#print(iso, introns)
				for r in bamObj.fetch(chrom, don - 20, acc + 20):
					op = check_read_compatible(bamObj, r, op, don, acc, introns, strand, args.l)
					
				##print(len(op))
				ncountAdj = 1.0
				for irn in op:
					opi = op[irn]
					num_reads = len(opi) 
					
					rnames = [r.query_name for r in opi]
					counts = Counter(rnames)
					
					nh = [r.get_tag("NH") for r in opi]
					nowcount = np.round(sum([1.0/(counts[rnames[k]]*nh[k]) for k in range(num_reads)]), 2)
					eventsFcounts[eid]["isos"][iso][irn].append(nowcount)

					ncountAdj *= max(nowcount, 1.0) 
					
				#print(ncountAdj, len(introns))
				ncountAdj = ncountAdj**(1/len(introns))
				events[eid]["iso"][iso][-1].append(np.round(ncountAdj, 2))
			

opf = os.path.join(args.o, args.i+".countsAdj.out")
with open(opf, "w") as oph:
	oph.write("gene_id\tfeature_id\t"+"\t".join(args.bams)+"\n")
	for eid in events:
		for iso in events[eid]["iso"]:
			oph.write(eid+"\t"+iso+"\t"+"\t".join([str(e) for e in events[eid]["iso"][iso][-1]])+"\n")


fopf = os.path.join(args.o, args.i+".feat.countsAdj.out")
with open(fopf, "w") as oph:
	oph.write("gene_id\tfeature_id\tjunctions\t"+"\t".join(args.bams)+"\n")
	for eid in eventsFcounts:
		for iso in eventsFcounts[eid]["isos"]:
			jns, jncounts = [], []
			
			for jn in eventsFcounts[eid]["isos"][iso]:
				jns.append(str(jn))

			for i in range(len(bams)):
				countstr = ""
				for jn in eventsFcounts[eid]["isos"][iso]:
					countstr += str(eventsFcounts[eid]["isos"][iso][jn][i]) + ","
				jncounts.append(countstr[:-1])

			oph.write(eid+"\t"+iso+"\t"+"_".join(jns)+"\t"+"\t".join(jncounts)+"\n")
