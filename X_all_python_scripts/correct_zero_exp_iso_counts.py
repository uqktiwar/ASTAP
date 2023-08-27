import sys
import os

"""
Feat file format:
   0			1							2				3					  4						5					6					7
isoform	HSC.C002UUB1_S001QBB1_S0025CB1	HSC.EXP	CLP.C002UUB3_S000SBB3_S001QBB3	CLP.EXP	CMP.C07002T2_C07015T3_C12001RP3	CMP.EXP	GMP.C002UUB4_C07002T3_C12001RP4	GMP.EXP	EB.C006MUB1_C006NSB2_C006POB1	EB.EXP	MEP.C000WYB3_C07002T4_C07015T4_C12001RP2	MEP.EXP	MK.C006NSB1_S0018AB4_S00198B6	MK.EXP	ANY.3.EXP
D2.event1.1	10.0_22.0_0.0	0	0.0_0.0_0.0	0	0.0_0.0_0.0	0	0.0_0.0_0.0	0	196.0_190.0_259.5	1	6.0_0.0_0.0_0.0	0	99.5_21.0_24.0	1	1
D2.event1.2	0.0,0.0_10.5,3.0_1.0,0.0	0	0.0,0.0_0.0,0.0_0.0,0.0	0	0.0,0.0_0.0,0.0_0.0,0.0	0	0.0,0.0_0.0,0.0_0.0,0.0	0	32.0,28.0_36.0,42.0_96.0,83.0	1	0.0,0.0_0.0,0.0_0.0,0.0_0.0,0.0	0	17.0,16.0_20.0,15.0_12.0,6.0	1	1
D2.event2.1	820.5_32.0_53.0	1	241.0_0.0_91.0	0	0.0_37.0_159.0	0	212.0_32.0_223.0	1	0.0_0.0_0.0	0	33.0_1.0_0.0_0.0	0	3.0_0.0_0.0	0	1


Iso File format:
gene_id	feature_id	C002UUB1	S001QBB1	S0025CB1	C002UUB3	S000SBB3	S001QBB3	C07002T2	C07015T3	C12001RP3	C002UUB4	C07002T3	C12001RP4	C006MUB1	C006NSB2	C006POB1	C000WYB3	C07002T4	C07015T4	C12001RP2	C006NSB1	S0018AB4	S00198B6
D2.event1	D2.event1.1	10.0	22.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	196.0	190.0	259.5	6.0	1.0	1.0	1.0	99.5	21.0	24.0
D2.event1	D2.event1.2	1.0	5.61	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	1.0	29.93	38.88	89.26	1.0	1.0	1.0	1.0	16.49	17.32	8.49
"""

def proc_feat_file(featfilepath):

	iso_all_feat_noexp = {}
	
	samples = []
	count_cols_idx = []
	with open(featfilepath, "r") as ffh:
		for line in ffh:
			if "isoform" in line:
				rec = line.strip().split("\t")
				i = 1
				for col in rec[1:-1]:
					if "EXP" in col:
						i += 1
						continue
										
					subsamples = col.split(".")[-1].split("_")
					samples.append(subsamples)
					count_cols_idx.append(i)
					i += 1
				
				continue
		
			rec = line.strip().split("\t")
			if rec[0] not in iso_all_feat_noexp:
				iso_all_feat_noexp[rec[0]] = {}
			
			# samples = [[C002UUB1, S001QBB1, S0025CB1], [C002UUB1, S001QBB1, S0025CB1], ...] 
			# count_cols_idx = [1, 3, 5, 7 ...]
			for (si, idx) in zip(samples, count_cols_idx): 
				fcounts = rec[idx].split("_")  ## ["0.0,0.0", "10.5,3.0", "1.0,0.0" ...]
				isoexp = [1 if fci.replace(",", "").replace("0", "").replace(".", "") == "" else 0 for fci in fcounts]

				for (sii, expi) in zip(si, isoexp):
					iso_all_feat_noexp[rec[0]][sii] = expi
			
	return iso_all_feat_noexp


			
iso_all_feat_noexp = proc_feat_file(featfilepath = os.path.join(sys.argv[1], sys.argv[2]))

isofilepath = os.path.join(sys.argv[1], sys.argv[3])
opfilepath = os.path.join(sys.argv[1], sys.argv[3]+".zero")

with open(opfilepath, "w") as oph:
	with open(isofilepath, "r") as inph:
		samples = []
		for line in inph:
			if "gene_id" in line:
				samples = line.strip().split("\t")[2:]
				oph.write(line)
				continue
			
			rec = line.strip().split("\t")

			fid = rec[1]
			newrec = rec.copy()

			for (i, s) in zip(range(2, len(rec)), samples):
				if iso_all_feat_noexp[fid][s]:
					newrec[i] = "0"
			
			oph.write("\t".join(newrec)+"\n")










