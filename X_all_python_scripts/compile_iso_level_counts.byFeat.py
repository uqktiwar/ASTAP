import numpy as np
import sys
import argparse as ap
import os

"""
###################### 
WHILE THIS WAS ORIGINALLY WRITTEN FOR PROCESSING BOTH NOTIR AND IR EVENTS FROM THE BLUEPRINT DATASET 
WE ONLY USE FOR IR EVENTS FOR THE NEURAL CREST 
#######################
"""

def map_iso_to_features(ffpath, etype, iso_to_feature = None):
	## this function retruns a dictionary of isoform
	## ids mapped to their features (only junctions for NOTIR events)
	## (junctions and retained introns for IR events)
	## The IR isoforms will contain all overlapping ret. introns
	## but the apt one for quantification will be filtered out after using the RI quant file 
	
	if iso_to_feature == None:
		iso_to_feature = {}

	if etype.upper() == "NOTIR":
		with open(ffpath, "r") as inph:
			for line in inph:
				rec = line.strip().split("\t")

				chrom, strand, iso = rec[0].replace("chr", ""), rec[5], rec[3]
				s, e = int(rec[1]), int(rec[2])
				
				estarts, elens = rec[11].split(",")[:-1], rec[10].split(",")[:-1]

				exons = [(s+int(es), s+int(es)+int(el)) for (es, el) in zip(estarts, elens)]

				numex = len(exons)

				introns = [(exons[i][1], exons[i+1][0]) for i in range(numex -1)]

				iso_to_feature[iso] = {"chrom":chrom, "strand":strand, "jns":introns, "exp":{}}
	
		return iso_to_feature

	if etype.upper() == "IR":
		with open(ffpath, "r") as inph:
			for line in inph:

				rec = line.split("\t")
				if rec[0] == "gene_id":
					continue
				
				chrom, strand, iso = rec[2].replace("chr", ""), rec[3], rec[1]
				jns, ri = [], []

				if rec[7] != "":
					jns = [tuple([int(ji) for ji in j.replace("(", "").replace(")", "").split(", ")]) for j in rec[7].split("),(")] ##(6648119, 6648256),(6648337, 6648502),(6648975, 6649340) >> ["(6648119, 6648256", "6648337, 6648502", "6648975, 6649340)"] >> [(6648119, 6648256), (6648119, 6648256) .. ]
				
				if rec[8].strip() != "":
					ri = list(set([tuple([int(ji) for ji in j.replace("(", "").replace(")", "").split(", ")]) for j in rec[8].split("),(")]))
				
				iso_to_feature[iso] = {"chrom":chrom, "strand":strand, "jns":jns, "ris":ri, "exp":{}}
			
		return iso_to_feature


def sample_info(pdpath):

	samples = {}
	
	with open(pdpath, "r") as pdh:
		for line in pdh:
			rec = line.strip().split("\t")
			
			if rec[0].upper() == "SAMPLE":
				continue
			
			samples[rec[0]] = rec[1:] ## map sample id to it's celltype (group) and sub-directory with count files
	
	return samples


def assemble_sampwise_jn_counts(masterdir, cfp, samples, jn_to_count = None):
	"""
	masterdir = master directory in which sub-directories containing countfiles are placed 
	cfp = count file prefix to filter files in the countfile directories
	samples = samples dictionary
	"""
	if jn_to_count == None:
		jn_to_count = {}
	
	countdirs = list(set([sv[-1] for sv in list(samples.values())]))
	
	for cd in countdirs:
		dpath = os.path.join(masterdir, cd)
		
		files = os.listdir(dpath)

		for f in files:
			if cfp in f:
				fpath = os.path.join(dpath, f)
				celltypes, samplekeys = None, None
				with open(fpath, "r") as inph:
					for line in inph:

						rec = line.strip().split("\t")

						if "strand" in rec:
							samplekeys = [r.replace(".count", "") for r in rec[5:]]
							celltypes = [samples[sk][0] for sk in samplekeys]
							continue
						
						jnidx = rec[4].replace("chr", "")

						if jnidx not in jn_to_count:
							jn_to_count[jnidx] = {}
						
						for (ct, sk, jc) in zip(celltypes, samplekeys, rec[5:]):
							if ct not in jn_to_count[jnidx]:
								jn_to_count[jnidx][ct] = {}
							
							if sk not in jn_to_count[jnidx][ct]:	
								jn_to_count[jnidx][ct][sk] = 0.0
						
							jn_to_count[jnidx][ct][sk] = float(jc)
						
	return jn_to_count


def assemble_sampwise_ri_counts(masterdir, cfp, samples, ri_to_count = None):

	"""
	 0	   1			2				3				4			5				6				7				8				9				10					11				12				13					14			15				16				17
    chr	strand	ret.intron.start	ret.intron.end	isoforms	unqiue.regions	C07002T4.count	C07002T4.quantl	C07002T4.stdev	C12001RP2.count	C12001RP2.quantl	C12001RP2.stdev	C07015T4.count	C07015T4.quantl	C07015T4.stdev	C000WYB3.count	C000WYB3.quantl	C000WYB3.stdev
	chr1	+	6648256	6648337	D3.event1.1	[6648256, 6648337]	0.0	(0.0, 0.0, 0.0)	0.7817847011845711	0.0	(0.0, 0.0, 0.0)	0.26764794307010864	0.0	(0.0, 0.0, 0.0)	0.22084621999997922	0.0	(0.0, 0.0, 0.0)	0.1551827788885992
	"""
	if ri_to_count == None:
		ri_to_count = {}
	
	countdirs = list(set([sv[-1] for sv in list(samples.values())]))
	
	for cd in countdirs:
		dpath = os.path.join(masterdir, cd)
		
		files = os.listdir(dpath)

		for f in files:
			if cfp in f:
				fpath = os.path.join(dpath, f)
				celltype, samplekeys, countcolz = None, None, None
				with open(fpath, "r") as inph:
					for line in inph:

						rec = line.strip().split("\t")

						if "strand" in rec:
							samplekeys = list(set([r.replace(".count", "") for r in rec[6:] if ".count" in r]))
							countcolz = [rec.index(sk+".count") for sk in samplekeys]# [ci for ci in range(6, len(rec)) if rec[ci].split(".")[-1] == "count"] ## get the indexes of the count columns in the ri quant file
							celltypes = [samples[sk][0] for sk in samplekeys]
							continue
						
						ridx = "__".join([rec[0].replace("chr", ""), rec[1], rec[2], rec[3]])   ## 1__+__783186__787306

						if ridx not in ri_to_count:
							ri_to_count[ridx] = {}
						
						for (ct, sk, jc) in zip(celltypes, samplekeys, countcolz):
							if ct not in ri_to_count[ridx]:
								ri_to_count[ridx][ct] = {}
							
							if sk not in ri_to_count[ridx][ct]:	
								ri_to_count[ridx][ct][sk] = 0.0
							
							if rec[5].upper() != "NO_UNQ_REGIONS":
								ri_to_count[ridx][ct][sk] = float(rec[jc])
						
	return ri_to_count

def remove_low_exp_iso_and_features(e2gpath, iso_to_features, etype):
	"""
	e2g format:
	 0	 1		2		3			4		  5			6		 7	  8		  9			10			11			12			13		14	15
	EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
	D3.event1	1	+	6648120	6649340	D3.event1.1,D3.event1.2,D3.event1.3	,6648256^6648338-6648502^6648816-6648904^6648976-,6648256^6648338-6648502^6648976-	MSTRG.388	ENSG00000204859	ZBTB48	6648502	6648976	D3.event1.2,D3.event1.3	6648816-6648904^,	3	2
	D2.event2	1	+	7837220	7838229	D2.event2.1,D2.event2.2	,7837378^7838178-	MSTRG.462	ENSG00000049245	VAMP3	7837220	7838229	D2.event2.1,D2.event2.2	,7837378^7838178-	2	2
	"""
	new_iso_to_features = {}
	with open(e2gpath, "r") as inph:
		for line in inph:

			rec = line.strip().split("\t")
			
			if rec[0] == "EID":
				continue
			
			isos_exp = rec[12].split(",")

			for iso in isos_exp:
				try:
					new_iso_to_features[iso] = iso_to_features[iso]
			
					if (rec[3], rec[4]) != (rec[10], rec[11]):
						nf0, nf1 = int(rec[10]), int(rec[11])
						
						rmjns = []
						for jn in new_iso_to_features[iso]["jns"]:
							if (jn[0] + 1) < nf0 or jn[1] <= nf0:
								rmjns.append(jn)
								
							if (jn[0] + 1) >= nf1 or jn[1] > nf1:
								rmjns.append(jn)
							
						new_iso_to_features[iso]["jns"] = [jn for jn in new_iso_to_features[iso]["jns"] if jn not in rmjns]
							
						if etype == "IR":
							rmris = []
							for ri in new_iso_to_features[iso]["ris"]:
								if (ri[0] + 1) < nf0 or ri[1] <= nf0:
									rmris.append(ri)
								
								if (ri[0] + 1) >= nf1 or ri[1] > nf1:
									rmris.append(ri)
								
							new_iso_to_features[iso]["ris"] = [ri for ri in new_iso_to_features[iso]["ris"] if ri not in rmris]
				except KeyError:
					break
	
	return new_iso_to_features


def add_expression_of_isoform_features(isoFeatures, samples, jnCounts, riCounts = None):
	celltypes = {}

	for sk in samples:
		ct = samples[sk][0]

		if ct not in celltypes:
			celltypes[ct] = []
		
		celltypes[ct].append(sk)
	

	for iso in isoFeatures:  ##iso_to_feature[iso] = {"chrom":chrom, "strand":strand, "jns":introns, "exp":{}}
		##print(isoFeatures[iso])
		for ct in celltypes:
			feature_counts_in_cell = []
			
			for sk in celltypes[ct]:
				
				feature_counts_in_cell_sample = []
				##print(isoFeatures[iso]["jns"])
				for jn in isoFeatures[iso]["jns"]:

					jnidx = "__".join([isoFeatures[iso]["chrom"], isoFeatures[iso]["strand"], str(jn[0]), str(jn[1])])
					
					feature_counts_in_cell_sample.append(max(jnCounts[jnidx][ct][sk], 1.0))
				
				if riCounts != None:
					for jn in isoFeatures[iso]["ris"]:

						jnidx = "__".join([isoFeatures[iso]["chrom"], isoFeatures[iso]["strand"], str(jn[0]), str(jn[1])])
						try:
							feature_counts_in_cell_sample.append(max(riCounts[jnidx][ct][sk], 1.0))
						except:
							pass
			
				gm_count = np.round(np.prod(feature_counts_in_cell_sample)**(1.0/len(feature_counts_in_cell_sample)), 2)
				feature_counts_in_cell.append(str(gm_count))

			isoFeatures[iso]["exp"][ct] = feature_counts_in_cell

	return (isoFeatures, celltypes)


parser = ap.ArgumentParser()
parser.add_argument("-w", "--w", help = "Path to working directory where input and output files are placed. Defaults to PWD.", default = os.getcwd())
parser.add_argument("-mjd1", "--masterDir1", help = "Full path to master directory containing sample level sub-directories with jn counts. Defaults to PWD. Must be same for JN counts and RI counts.", default = os.getcwd())
parser.add_argument("-mjd2", "--masterDir2", help = "Full path to master directory containing sample level sub-directories with jn counts. Defaults to PWD. Must be same for JN counts and RI counts.", default = os.getcwd())
parser.add_argument("-pdf1", "--pdataFile1", help = "Name of the pdata file containing sample details. Must be placed in working directory.")
parser.add_argument("-pdf2", "--pdataFile2", help = "Name of the pdata file containing sample details. Must be placed in working directory. Default = None", default = None)
parser.add_argument("-cfpJn1", "--jnCountFilePrefix1", help = "Common prefix of the file name within the sub-drectories that contains the junction counts")
parser.add_argument("-cfpJn2", "--jnCountFilePrefix2", help = "Common prefix of the file name within the sub-drectories that contains the junction counts. Default = None", default = None)
parser.add_argument("-cfpRi", "--riCountFilePrefix", help = "Common prefix of the file name within the sub-drectories that contains the retained intron counts. Optional", nargs = "?")
parser.add_argument("-iff1", "--inputFeatureFile1", help = "Name of input file containing isoform level bed records or features for the events. File must be in the working directory.")
parser.add_argument("-iff2", "--inputFeatureFile2", help = "Name of input file containing isoform level bed records or features for the events. File must be in the working directory. Default = None", default = None)
parser.add_argument("-e2gf", "--e2gFile", help = "Name of the event2gene file containing info about expressed isoforms and new flanks for each event.")
parser.add_argument("-etype", "--etype", help = "Type of event. IR or NOTIR. Case insensitive. Default = NOTIR", default = "notir")
parser.add_argument("-op", "--outPrefix", help = "Prefix for output file name.")
args = parser.parse_args()					

ffpath1 = os.path.join(args.w, args.inputFeatureFile1)
iso_to_feature = map_iso_to_features(ffpath1, args.etype)

pdpath1 = args.pdataFile1
samples1 = sample_info(pdpath1)

masterdir1, masterdir2 = os.path.join(args.w, args.masterDir1), os.path.join(args.w, args.masterDir2)
jn_to_count = assemble_sampwise_jn_counts(masterdir1, args.jnCountFilePrefix1, samples1)
ri_to_count = None

if args.etype.upper() == "IR":
	ri_to_count = assemble_sampwise_ri_counts(masterdir1, args.riCountFilePrefix, samples1)

if args.inputFeatureFile2 != None:
	ffpath2 = os.path.join(args.w, args.inputFeatureFile2)
	iso_to_feature = map_iso_to_features(ffpath2, args.etype, iso_to_feature)

	pdpath2 = args.pdataFile2
	samples2 = sample_info(pdpath2)

	jn_to_count = assemble_sampwise_jn_counts(masterdir2, args.jnCountFilePrefix2, samples2, jn_to_count)

	if args.etype.upper() == "IR":
		ri_to_count = assemble_sampwise_ri_counts(masterdir2, args.riCountFilePrefix, samples2, ri_to_count)

"""
Feature Count dictionaries structure
1__+__6648256__6648337
{
 'EB': {'C006MUB1': 82.0, 'C006NSB2': 39.0, 'C006POB1': 68.0}, 
 'MEP': {'C07002T4': 183.0, 'C12001RP2': 32.83, 'C07015T4': 85.0, 'C000WYB3': 95.5}, 
 'CLP': {'S001QBB3': 16.0, 'S000SBB3': 0.0, 'C002UUB3': 187.0}, 
 'GMP': {'C002UUB4': 209.0, 'C07002T3': 269.0, 'C12001RP4': 278.0}, 
 'MPP': {'C07015T1': 575.0, 'S001FXB3': 0.0}, 
 'MK': {'C006NSB1': 127.0, 'S0018AB4': 85.0, 'S00198B6': 47.0}, 
 'HSC': {'S001QBB1': 29.0, 'S0025CB1': 10.0, 'C002UUB1': 102.0}, 
 'CMP': {'C07002T2': 145.33, 'C12001RP3': 122.83, 'C07015T3': 229.17}
}
"""
e2gpath = os.path.join(args.w, args.e2gFile)
if args.inputFeatureFile2 != None:
	iso_to_feature = remove_low_exp_iso_and_features(e2gpath, iso_to_feature, args.etype.upper())

iso_to_feature, celltypes = add_expression_of_isoform_features(iso_to_feature, samples1, jn_to_count, ri_to_count)
countsOpath = os.path.join(args.w, args.outPrefix + ".IsoExprsn.out")

with open(countsOpath, "w") as coph:
	headr = ["gene_id", "feature_id"]
	
	for ct in celltypes:
		for sk in celltypes[ct]:
			headr.append(sk)
	
	coph.write("\t".join(headr) + "\n")

	for iso in iso_to_feature:
		eid = ".".join(iso.split(".")[:-1])
		orec = [eid, iso]

		for ct in celltypes:
			orec.extend(iso_to_feature[iso]["exp"][ct])
		
		coph.write("\t".join(orec) + "\n")

