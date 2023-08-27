import pysam as pys
import sys
import numpy as np
import os

"""
1	1148473	1149042	D2.event739.1	0	-
1	1670692	1670929	D2.ovevent113.1	0	-
1	1683836	1683880	D2.event760.1	0	-
1	1718110	1718809	D4.event68.1	0	-
1	1718492	1718760	D4.event68.1	0	-
1	1718492	1718764	D4.event68.1	0	-
1	26800076	26800193	D5.event4.5	0	+
1	26800076	26800193	D6.event1.1	0	+
"""

def pre_proc_ret_jns(ripath):
	"""
	Use this function remove redundant introns 
	"""
	retJns = {}
	
	with open(ripath, "r") as rih:
		for line in rih:
			rec = line.strip().split("\t")
			
			idx = "__".join([rec[0], rec[-1], rec[1], rec[2]]) ## chr, strand, start (0-b), end (1-b)

			if idx not in retJns:
				retJns[idx] = []
			
			retJns[idx].append(rec[3])
	
	return retJns


retJns = pre_proc_ret_jns(os.path.join(sys.argv[1], sys.argv[2]))
opfile = os.path.join(sys.argv[1], sys.argv[3])
check_strand = int(sys.argv[4])

tbx = pys.TabixFile(os.path.join(sys.argv[1], "asg.exon_database.bed.gz"))

with open(opfile, "w") as outh:
	outh.write("\t".join(["chr", "strand", "ret.intron.start", "ret.intron.end", "tx.ids", "unqiue.regions", "ovlap.exons"]) + "\n")
	for ri in retJns:
		chrom, strand, ris, rie = ri.split("__")
		ris, rie = int(ris), int(rie)
		
		ricoords = list(range(ris, rie))
		ribases = np.zeros(len(ricoords))
		ovlpeid = []
		for row in tbx.fetch(chrom, ris, rie):
			##
			rec = str(row).split("\t") ## [0'1', 1'1148371', 2'1149165', 3'exon.186', 40', 5'-']

			if check_strand and rec[-1] != strand:
				continue

			ovlap_coords = [max(int(rec[1]), ris), min(int(rec[2]), rie)]
			if ovlap_coords == [ris, rie]: ## dont consider the exons formed by retaining the intron, which would contain the full intron 
				continue

			a, b = ricoords.index(ovlap_coords[0]), ricoords.index(ovlap_coords[1] - 1) ## both are zero based now
			ribases[a:b + 1] = 1
			ovlpeid.append(rec[3])
		
		ribases = list(ribases)
		non_overlapping_regions = []
		
		interval = []
		for i in range(len(ribases)):
			if ribases[i] == 0:
				if interval == []:
					interval = [ricoords[i], 0]
				else:
					interval[1] = ricoords[i] + 1  ## 0, 1 based 
			
			if ribases[i] == 1 and interval != []:
				non_overlapping_regions.append(interval)
				interval = []
		
		if interval != [] and (interval not in non_overlapping_regions):
			if interval[1] == 0:
				interval[1] = ricoords[-1] + 1
			non_overlapping_regions.append(interval)
		
		if non_overlapping_regions == []:
			outh.write("\t".join([str(e) for e in [chrom, strand, ris, rie, ",".join(retJns[ri]), "No_Unq_Regions"]]) + "\n")
			continue

		non_overlapping_regions = "_".join([str(e) for e in non_overlapping_regions])
		outh.write("\t".join([str(e) for e in [chrom, strand, ris, rie, ",".join(retJns[ri]), non_overlapping_regions, ",".join(ovlpeid)]]) + "\n")  ### 1111111111111111111111111110
		
		






				


		



			



