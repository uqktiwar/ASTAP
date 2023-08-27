import sys
import argparse as ap
import os


"""
HMP IR feature file format
   0		1		  2		  3			4			5	  6			7				8
gene_id	feature_id	chrom	strand	ev.start	ev.end	exons	junctions	retained_introns
D3.event1	D3.event1.1	1	+	6648119	6649340	(6648119, 6649340)		(6648256, 6648337),(6648502, 6648815),(6648904, 6648975),(6648256, 6648337),(6648502, 6648975)
D3.event1	D3.event1.2	1	+	6648119	6649340	(6648119, 6648256),(6648337, 6648502),(6648815, 6648904),(6648975, 6649340)	(6648256, 6648337),(6648502, 6648815),(6648904, 6648975)	
D3.event1	D3.event1.3	1	+	6648119	6649340	(6648119, 6648256),(6648337, 6648502),(6648975, 6649340)	(6648256, 6648337),(6648502, 6648975)	


HMP NotIR BED format
0		1		2		3		4	5		6		7	8	9	  10		11
1	783136	787356	D2.event1.1	1	+	783136	787356	1	2	50,50,	0,4170,
1	783136	787356	D2.event1.2	1	+	783136	787356	1	3	50,119,50,	0,1727,4170,
1	764434	787356	D3.event2.1	1	+	764434	787356	1	2	50,50,	0,22872,
1	764434	787356	D3.event2.2	1	+	764434	787356	1	3	50,153,50,	0,18599,22872,

JN Counts file -- For NOTIR events
   0		1			2		  3		   4	  5		  6
gene_id	feature_id	junctions	24701	24703	24704	24705
D2.event1	D2.event1.1	(321264, 322037)	3.08	1.33	1.92	2.42
D2.event1	D2.event1.2	(321290, 322037)	1.0	0	0	1.0
D2.event2	D2.event2.1	(783186, 787306)	10.0	12.0	5.0	7.0

JN Count file format -- Only for IR events
 0		1			2		3	4			5			6				7
chr	jn.start	jn.end	strand	id	C006MUB1.count	C006NSB2.count	C006POB1.count
1	+	783186	787306	1__+__783186__787306	196.0	190.0	259.5

RI Count file format
 0		1			2				3				4			5				6				7				8				9				10					11				12				13					14			15				16				17
chr	strand	ret.intron.start	ret.intron.end	isoforms	unqiue.regions	C07002T4.count	C07002T4.quantl	C07002T4.stdev	C12001RP2.count	C12001RP2.quantl	C12001RP2.stdev	C07015T4.count	C07015T4.quantl	C07015T4.stdev	C000WYB3.count	C000WYB3.quantl	C000WYB3.stdev
chr1	+	6648256	6648337	D3.event1.1	[6648256, 6648337]	0.0	(0.0, 0.0, 0.0)	0.7817847011845711	0.0	(0.0, 0.0, 0.0)	0.26764794307010864	0.0	(0.0, 0.0, 0.0)	0.22084621999997922	0.0	(0.0, 0.0, 0.0)	0.1551827788885992
chr1	+	6648502	6648975	D3.event1.1	[6648502, 6648746]_[6648904, 6648975]	0.0	(0.0, 0.0, 0.0)	1.9572503478508463	0.0	(0.0, 0.0, 0.0)	0.41492075636679826	0.0	(0.0, 0.0, 0.0)	0.9175042555524832	0.0	(0.0, 0.0, 0.0)	0.550123099442467

pdata file format
  0			1
Sample	CellType
S001QBB1	HSC
S0025CB1	HSC
C002UUB1	HSC
"""


def map_iso_to_features(fcpath, etype):
	## this function retruns a dictionary of isoform
	## ids mapped to their features (only junctions for NOTIR events)
	## (junctions and retained introns for IR events)
	## The IR isoforms will contain all overlapping ret. introns
	## but the apt one for quantification will be filtered out after using the RI quant file 

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


def assemble_sampwise_jn_counts(masterdir, cfp, samples, etype, isofeatures = None):
	"""
	masterdir = master directory in which sub-directories containing countfiles are placed 
	cfp = count file prefix to filter files in the countfile directories
	samples = samples dictionary
	"""
	jn_to_count = {}
	countdirs = list(set([sv[-1] for sv in list(samples.values())]))
	if etype == "NOTIR":
		"""
		   0		1			2		  3		  4		  5		  6
		gene_id	feature_id	junctions	24701	24703	24704	24705
		D2.event1	D2.event1.1	(321264, 322037)	3.08	1.33	1.92	2.42
		D2.event1	D2.event1.2	(1244100, 1244294)_(1244352, 1244465)	1.0	0	0	1.0
		"""
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

							if "gene_id" in rec:
								samplekeys = [r.replace(".count", "") for r in rec[3:]]
								celltypes = [samples[sk][0] for sk in samplekeys]
								continue
							
							chrom, strand = isofeatures[rec[1]]["chrom"], isofeatures[rec[1]]["strand"]
							
							for jn in rec[2].split("_"):
								jnidx = "__".join([chrom, strand] + jn.replace("(", "").replace(")", "").split(", ")) ## chrom__strand__jn0__jn1

								if jnidx not in jn_to_count:
									jn_to_count[jnidx] = {}
								
								for (ct, sk, jc) in zip(celltypes, samplekeys, rec[3:]):
									if ct not in jn_to_count[jnidx]:
										jn_to_count[jnidx][ct] = {}
									
									if sk not in jn_to_count[jnidx][ct]:	
										jn_to_count[jnidx][ct][sk] = 0.0
								
									jn_to_count[jnidx][ct][sk] = float(jc.split(",")[rec[2].split("_").index(jn)])
		return jn_to_count	
		

	if etype == "IR":
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


def assemble_sampwise_ri_counts(masterdir, cfp, samples):

	"""
	 0	   1			2				3				4			5				6				7				8				9				10					11				12				13					14			15				16				17
    chr	strand	ret.intron.start	ret.intron.end	isoforms	unqiue.regions	C07002T4.count	C07002T4.quantl	C07002T4.stdev	C12001RP2.count	C12001RP2.quantl	C12001RP2.stdev	C07015T4.count	C07015T4.quantl	C07015T4.stdev	C000WYB3.count	C000WYB3.quantl	C000WYB3.stdev
	chr1	+	6648256	6648337	D3.event1.1	[6648256, 6648337]	0.0	(0.0, 0.0, 0.0)	0.7817847011845711	0.0	(0.0, 0.0, 0.0)	0.26764794307010864	0.0	(0.0, 0.0, 0.0)	0.22084621999997922	0.0	(0.0, 0.0, 0.0)	0.1551827788885992
	"""
	
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
							samplekeys = [r.replace(".count", "") for r in rec[6:] if "count" in r]
							countcolz = [ci for ci in range(6, len(rec)) if rec[ci].split(".")[-1] == "count"] ## get the indexes of the count columns in the ri quant file
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


def add_expression_of_isoform_features(isoFeatures, samples, jnCounts, riCounts = None):
	celltypes = {}

	for sk in samples:
		ct = samples[sk][0]

		if ct not in celltypes:
			celltypes[ct] = []
		
		celltypes[ct].append(sk)
	
	print(celltypes)
	for iso in isoFeatures:  ##iso_to_feature[iso] = {"chrom":chrom, "strand":strand, "jns":introns, "exp":{}}

		for ct in celltypes:
			feature_counts_in_cell = []
			
			for sk in celltypes[ct]:
				
				feature_counts_in_cell_sample = []
	
				for jn in isoFeatures[iso]["jns"]:

					jnidx = "__".join([isoFeatures[iso]["chrom"], isoFeatures[iso]["strand"], str(jn[0]), str(jn[1])])
					
					feature_counts_in_cell_sample.append(jnCounts[jnidx][ct][sk])
				
				if riCounts != None:
					for jn in isoFeatures[iso]["ris"]:

						jnidx = "__".join([isoFeatures[iso]["chrom"], isoFeatures[iso]["strand"], str(jn[0]), str(jn[1])])
						try:
							feature_counts_in_cell_sample.append(riCounts[jnidx][ct][sk])
						except:
							pass
			
				feature_counts_in_cell.append([feature_counts_in_cell_sample, min(feature_counts_in_cell_sample)])

			isoFeatures[iso]["exp"][ct] = feature_counts_in_cell

	return (isoFeatures, celltypes)


parser = ap.ArgumentParser()
parser.add_argument("-w", "--w", help = "Path to working directory where input and output files are placed. Defaults to PWD.", default = os.getcwd())
parser.add_argument("-mjd", "--masterDir", help = "Full path to master directory containing sample level sub-directories with jn counts. Defaults to PWD. Must be same for JN counts and RI counts.", default = os.getcwd())
parser.add_argument("-pdf", "--pdataFile", help = "Name of the pdata file containing sample details. Must be placed in working directory.")
parser.add_argument("-cfpJn", "--jnCountFilePrefix", help = "Common prefix of the file name within the sub-drectories that contains the junction counts")
parser.add_argument("-cfpRi", "--riCountFilePrefix", help = "Common prefix of the file name within the sub-drectories that contains the retained intron counts. Optional", nargs = "?")
parser.add_argument("-iff", "--inputFeatureFile", help = "Name of or path within working dir to input file containing isoform level bed records or features for the events. File must be somewhere in the working directory.")
parser.add_argument("-ft", "--featureThresh", help = "Count threshold for features to be considered expressed. Default = 2", type = int, default = 2)
parser.add_argument("-mins", "--minSamples", help = "Number of samples to check for expression. Default = 3", type = int, default = 3)
parser.add_argument("-etype", "--etype", help = "Type of event. IR or NOTIR. Case insensitive. Default = NOTIR", default = "notir")
parser.add_argument("-op", "--outPrefix", help = "Prefix for output file name.")
args = parser.parse_args()					

ffpath = os.path.join(args.w, args.inputFeatureFile)
iso_to_feature = map_iso_to_features(ffpath, args.etype)
	
pdpath = os.path.join(args.w, args.pdataFile)
samples = sample_info(pdpath)
print(samples)
masterdir = args.masterDir
jn_to_count, ri_to_count = None, None

if args.etype.upper() == "NOTIR": 
	jn_to_count = assemble_sampwise_jn_counts(masterdir, args.jnCountFilePrefix, samples, args.etype.upper(), iso_to_feature)				

if args.etype.upper() == "IR":
	jn_to_count = assemble_sampwise_jn_counts(masterdir, args.jnCountFilePrefix, samples, args.etype.upper())
	ri_to_count = assemble_sampwise_ri_counts(masterdir, args.riCountFilePrefix, samples)

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
iso_to_feature, celltypes = add_expression_of_isoform_features(iso_to_feature, samples, jn_to_count, ri_to_count)

"""
isos = []
with open("isos", "r") as isoh:
	for line in isoh:
		isos.append(line.strip())
"""

countsOpath = os.path.join(args.w, args.outPrefix + ".FeatExprsn.out")

with open(countsOpath, "w") as coph:
	headr = ["isoform"]
	
	for ct in celltypes:
		headr.append(ct + "." + "_".join(celltypes[ct]))
		headr.append(ct + ".EXP")
	
	headr.append("ANY." + str(args.minSamples) + ".EXP")
	coph.write("\t".join(headr) + "\n")

	for iso in iso_to_feature:
		orec = [iso]
		master_nsamp = 0
		for ct in celltypes:
			
			nsamp = sum([1 if elem[-1] >= args.featureThresh else 0 for elem in iso_to_feature[iso]["exp"][ct]])
			countstring = "_".join([",".join([str(ec) for ec in elem[0]]) for elem in iso_to_feature[iso]["exp"][ct]])

			if nsamp >= args.minSamples:
				orec.extend([countstring, "1"])
			
			if nsamp < args.minSamples:
				orec.extend([countstring, "0"])

			master_nsamp += nsamp 
			
		if master_nsamp >= args.minSamples:
			orec.append("1")
		
		if master_nsamp < args.minSamples:
			orec.append("0")
		
		coph.write("\t".join(orec) + "\n")









