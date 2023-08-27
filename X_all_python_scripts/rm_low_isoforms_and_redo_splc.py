import sys
import re
from collections import Counter
"""
Arg1 - E2G FILE

Arg2 - Feat Exp file 
isoform	HSC.S001QBB1_S0025CB1_C002UUB1	HSC.EXP	MPP.S001FXB3_C07015T1	MPP.EXP	CMP.C12001RP3_C07015T3_C07002T2	CMP.EXP	CLP.C002UUB3_S000SBB3_S001QBB3	CLP.EXP	GMP.C002UUB4_C07002T3_C12001RP4	GMP.EXP	MEP.C000WYB3_C12001RP2_C07002T4_C07015T4	MEP.EXP	MK.C006NSB1_S0018AB4_S00198B6	MK.EXP	EB.C006MUB1_C006NSB2_C006POB1	EB.EXP
D2.event1.1	22.0_0.0_10.0	1	0.0_0.0	0	0.0_0.0_0.0	0	0.0_0.0_0.0	0	0.0_0.0_0.0	0	6.0_0.0_0.0_0.0	0	99.5_21.0_26.0	1	196.0_190.0_259.5	1

Arg3 - Op File name
"""

def compile_exp_info(featExpPath):
	
	iso_exp = {}
	exp_colz = []
	with open(featExpPath, "r") as inph:
		for line in inph:

			rec = line.strip().split("\t")

			if "isoform" in rec:
				##print(range(len(rec)))
				exp_colz = [i for i in range(len(rec)) if rec[i].split(".")[-1] == "EXP"]
				continue
			
			eid = ".".join(rec[0].split(".")[:-1])

			if eid not in iso_exp:
				iso_exp[eid] = {}

			iso_exp[eid][rec[0]] = 0
			exp_vec = [int(rec[i]) for i in exp_colz]
			
			if sum(exp_vec):
				iso_exp[eid][rec[0]] = 1
			
	return iso_exp

def reject_for_flanks(s, newsplc):
	s = str(s)
	for splc in newsplc:
		if s not in splc:
			return False
		
	return True

def create_new_splchain(orig_splice_chain, strand, flanks, expSubD):
	"""
	D6.event21646	X	-	71694579	71791907
	71788734-71788604^71787880-71787739^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71700300-71700300^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710775^71709113-71709112^71708885-71708783^71700300-71700300^,71788734-71788604^71787880-71787739^71715118-71715006^71708891-71708783^71700300-71700300^,71715118-71715006^71708891-71708783^	
	MSTRG.18776	
	ENSG00000147099	
	HDAC8	71700300	71791907	
	D6.event21646.1,D6.event21646.2,D6.event21646.3,D6.event21646.4,D6.event21646.5,D6.event21646.6	
	71788734-71788604^71787880-71787739^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71700300-71700300^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710775^71709113-71709112^71708885-71708783^71700300-71700300^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71708891-71708783^71700300-71700300^71694579f,71715118-71715006^71708891-71708783^71694579f	6	6
	"""
	orig_splice_chain = orig_splice_chain.split(",") ## get the original splice chain for the event and split it into a list 
	new_splice_chain = []
	
	isos = list(expSubD.keys())
	isosno = [int(iso.split(".")[-1]) for iso in isos]
	
	isos_sorted = [x for _,x in sorted(zip(isosno,isos))]
	isos = isos_sorted

	isos_exp = [iso for iso in isos if expSubD[iso]]
	
	rd = len(isos_exp)

	if rd < 2:
		return None

	if strand == "+":
		new_splice_chain = [str(flanks[0]) + "f" + orig_splice_chain[i] + str(flanks[1]) + "f" for i in range(len(isos)) if expSubD[isos[i]]] ## create a new_splice_chain consiisting of the original splice chains 
																														   ## of the isoforms that are expressed acc to the read count threshold
	if strand == "-":
		new_splice_chain = [str(flanks[1]) + "f" + orig_splice_chain[i] + str(flanks[0]) + "f" for i in range(len(isos)) if expSubD[isos[i]]]
		##print(orig_splice_chain, strand, flanks, expSubD)
		
	
	## Create a list of all internal sites from the splice chains of the expressed transcripts (isoforms)
	all_sites = []  
	for iso in new_splice_chain:
		all_sites.extend([int(e) for e in re.split(r'[f\^\-\[\]]', iso)[:-1]])
	
	## Convert the list into a counter object to count occurence of each site and compare it with real dimension of the event 
	all_sites = list(set(all_sites))
	
	#print(isos)
	#print(all_sites)


	nf0, nf1 = sorted(all_sites)[0], sorted(all_sites)[-1]
	for s in sorted(all_sites):
		if not reject_for_flanks(s, new_splice_chain):
			break
		nf0 = s

	for s in sorted(all_sites)[::-1]:
		if not reject_for_flanks(s, new_splice_chain):
			break
		nf1 = s  ## at the start of the loop, nf0 = 0 all common sites will get added to cand[0] >> when the first and any subsequent variable site is encountered nf0 = 1 (max(0,1) or max(1,1)) >> when common sites are encountered after this they get added to cand[1]  

	nf0, nf1 = str(nf0), str(nf1)
	
	## For each internal site which occurs in all expressed isoforms set the smallest magnitude site as the new flank_0 and the max magnitude site as the new flank_1
	if strand == "+":
		try:
			new_splice_chain = [e[e.index(nf0) + len(nf0) + 1:] for e in new_splice_chain]  ## truncate the splice chains to start after the end of the 
		except:
			pass
		
		try:
			new_splice_chain = [e[:e.index(nf1)] for e in new_splice_chain] ## reduce all elements by the common new end flank  
		except:
			pass
	
	if strand == "-":  ## since the order of the coordinates in the splice chains are decreasing  37^32-28^ 
		##print(all_sites)
		##print(nf0, nf1)
		##print(new_splice_chain)
		try:
			new_splice_chain = [e[:e.index(nf0)] for e in new_splice_chain]
		except:
			pass
		
		try:
			
			new_splice_chain = [e[e.index(nf1) + len(nf1) + 1:] for e in new_splice_chain] ## reduce all elements by the common new end flank  
		except:
			pass
	
	new_splice_chain = ",".join(new_splice_chain)
	
	return (isos, isos_exp, nf0, nf1, new_splice_chain, rd)

e2gpath = sys.argv[1]
iso_exp = compile_exp_info(featExpPath = sys.argv[2])

with open(sys.argv[3], "w") as oph:
	oph.write("EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD\n")
	
	with open(e2gpath, "r") as e2h:
		for line in e2h:

			"""
			 0	 1	   2		3			4			5		 6	  7		 8
			EID	CHR	STRAND	FLANK.5P	FLANK.3P	SPLCHAIN	GID	EGID	HGNC
			event1	1	+	783186	787307	,784864-784982^	MSTRG.93	ENSG00000228794,ENSG00000234711,ENSG00000269308	LINC01128,TUBB8P11,AL645608.2
			"""
			
			rec = line.strip().split("\t")

			if rec[0] == "EID":
				continue
			
			dim = len(rec[5].split(","))

			eid = "D" + str(dim) + "." + rec[0]
			
			## oph.write("EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD\n")

			if dim < 3:
				orec = [eid] + rec[1:5] + [",".join(list(sorted(iso_exp[eid].keys())))] + rec[5:] + rec[3:5]  + [",".join(list(sorted(iso_exp[eid].keys()))), rec[5], str(dim), str(dim)]
				oph.write("\t".join(orec)+"\n")
				continue

			##print(eid)	
			outputt = create_new_splchain(orig_splice_chain = rec[5], strand = rec[2], flanks = [rec[3], rec[4]], expSubD = iso_exp[eid])
			if outputt == None:
				continue
			
			isos, isos_exp, nf0, nf1, new_splice_chain, rd = outputt
			orec = [eid] + rec[1:5] + [",".join(isos)] + rec[5:] + [nf0, nf1] + [",".join(isos_exp), new_splice_chain, str(dim), str(rd)]
			oph.write("\t".join(orec)+"\n")
			









