import numpy as np
import sys
import argparse as ap
import os

"""
   0		1		  2		  3		  4		  5
gene_id	feature_id	24701	24703	24704	24705
D2.event1	D2.event1.1	3.08	1.33	1.92	2.42
D2.event1	D2.event1.2	1.0	1.0	1.0	1.0
D2.event2	D2.event2.1	10.0	12.0	5.0	7.0
D2.event2	D2.event2.2	2.83	2.0	1.0	1.0
D3.event3	D3.event3.1	2.0	2.0	1.0	2.0
"""

def assemble_sampwise_jn_counts(masterdir, cfp, samples, iso_to_count = None, dusra = False):
	"""
	masterdir = master directory in which sub-directories containing countfiles are placed 
	cfp = count file prefix to filter files in the countfile directories
	samples = samples dictionary

	gene_id	feature_id	S01	S02	S03
	D2.event1	D2.event1.1	1099.0	1096.0	1136.0
	D2.event1	D2.event1.2	1160.17	1181.2	1207.19
	"""
	if iso_to_count == None:
		iso_to_count = {}
	
	allsamps = []
	for cd in samples:
		print(cd)
		dpath = os.path.join(masterdir, cd)
		
		try:
			files = os.listdir(dpath)
		except FileNotFoundError:
			continue

		samps = []
		for f in files:
			print(f)
			if cfp in f:
				fpath = os.path.join(dpath, f)
				#print(fpath)

				with open(fpath, "r") as inph:
					for line in inph:

						rec = line.strip().split("\t")
						if rec[0] == "gene_id":
							samps = rec[2:]
							allsamps.extend(samps)
							#if dusra:
							#	print(samps, allsamps)
							continue
						
						if rec[1] not in iso_to_count:
							iso_to_count[rec[1]] = {"eid":rec[0]}
						
						for (s,c) in zip(samps, rec[2:]):
							iso_to_count[rec[1]][s] = float(c)
							#if dusra:
							#	print(iso_to_count[rec[1]])

	
	for iso in iso_to_count:
		for s in allsamps:
			if s not in iso_to_count[iso]:
				iso_to_count[iso][s] = 0.0

						
	return (allsamps, iso_to_count)


parser = ap.ArgumentParser()
parser.add_argument("-w", "--w", help = "Path to working directory where input and output files are placed. Defaults to PWD.", default = os.getcwd())
parser.add_argument("-mjd1", "--masterDir1", help = "Full path to master directory containing subdirectories with isoform level counts for NotIR events.", default = os.getcwd())
parser.add_argument("-mjd2", "--masterDir2", help = "Full path to master directory containing subdirectories with requant NotIR event count files. Default = None", default = None)
parser.add_argument("-e2gf", "--e2gfile", help = "Name of final e2g file for filtering out events.")
parser.add_argument("-csubd", "--countSubDirs", help = "Names of sub-directories within master directories which contain iso level count files. Supply at least 1", nargs = "+")
parser.add_argument("-cfp1", "--jnCountFilePrefix1", help = "Common prefix of the file name within the sub-drectories in master directory 1")
parser.add_argument("-cfp2", "--jnCountFilePrefix2", help = "Common prefix of the file name within the sub-drectories in master directory 2. Default = None", default = None)
parser.add_argument("-op", "--outPrefix", help = "Prefix for output file name.")
args = parser.parse_args()					

masterdir1 = os.path.join(args.w, args.masterDir1)
allsamps, iso_to_count = assemble_sampwise_jn_counts(masterdir1, args.jnCountFilePrefix1, args.countSubDirs)				

if args.masterDir2 != None:
	masterdir2 = os.path.join(args.w, args.masterDir2)
	allsamps, iso_to_count = assemble_sampwise_jn_counts(masterdir2, args.jnCountFilePrefix2, args.countSubDirs, iso_to_count = iso_to_count, dusra = True)

countsOpath = os.path.join(args.w, args.outPrefix + ".IsoExprsn.out")

events2keep = []
with open(os.path.join(args.w, args.e2gfile), "r") as e2gh:
	for line in e2gh:
		if "EID" in line:
			continue
		if args.masterDir2 != None:
			events2keep.extend(line.strip().split("\t")[12].split(","))
		
		if args.masterDir2 == None:
			events2keep.extend(line.strip().split("\t")[5].split(","))
		"""
		 0	 1		2		3			4		 5			6		 7	  8		 9			10			11			12			13		14	15
		EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
		D2.event1	1	+	321032	322038	D2.event1.1,D2.event1.2	321264^,321290^	MSTRG.17	ENSG00000237094,ENSG00000250575,ENSG00000224813	RP4-669L17.10,RP4-669L17.4,RP4-669L17.8	321032	322038	D2.event1.1,D2.event1.2	321264^,321290^	2	2
		D2.event2	1	+	783186	787307	D2.event2.1,D2.event2.2	,784864-784982^	MSTRG.52	ENSG00000228794	LINC01128	783186	787307	D2.event2.1,D2.event2.2	,784864-784982^	2	2
		"""
 

with open(countsOpath, "w") as coph:
	headr = ["gene_id", "feature_id"] + allsamps
	coph.write("\t".join(headr) + "\n")

	for iso in iso_to_count:
		eid = iso_to_count[iso]["eid"]
		if iso not in events2keep:
			continue

		orec = [eid, iso]
		for s in allsamps:
			orec.append(str(iso_to_count[iso][s]))

		coph.write("\t".join(orec) + "\n")

