import sys

"""
Arg 1
Sample input Stringtie
chr1	StringTie	transcript	14363	29693	1000	-	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.3"; 
 0			1			2		  3		  4		  5		6	7	  8_0		8_1		 8_2			8_3			
chr1	StringTie	transcript	14363	29343	1000	-	.	gene_id "MSTRG.2"; transcript_id "MSTRG.2.4"; 

chr1	StringTie	transcript	14363	29370	1000	-	.	gene_id "MSTRG.2"; transcript_id "ENST00000423562"; gene_name "WASH7P"; ref_gene_id "ENSG00000227232"; 
 0			1			2		  3		  4		  5		6	7	  8_0		8_1		 8_2			8_3					8_4		8_5			8_6			8_7
chr1	StringTie	transcript	14363	24886	1000	-	.	gene_id "MSTRG.2"; transcript_id "ENST00000541675"; gene_name "WASH7P"; ref_gene_id "ENSG00000227232"; 

Sample Input GENCODE
 0         1      2					3	  4		5		6		7			8_0		8_1						8_2				8_3				  8_4				8_5								8_6		 8_7		8_8					8_9
chr1    HAVANA  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972.5_4"; transcript_id "ENST00000456328.2_1"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level 1; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2_4"; havana_transcript "OTTHUMT00000362751.1_1"; remap_num_mappings 1; remap_status "full_contig"; remap_target_status "overlap";



Sample Input Reference Ensembl
0		1						2		   3	   4	5	6	7	   8_0		8_1					8_2				8_3				8_4		8_5			8_6				8_7				8_8			8_9
1	processed_transcript	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
1	transcribed_unprocessed_pseudogene	transcript	11872	14412	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000515242"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-201"; transcript_source "ensembl";
1	transcribed_unprocessed_pseudogene	transcript	11874	14409	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000518655"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-202"; transcript_source "ensembl";
1	transcribed_unprocessed_pseudogene	transcript	12010	13670	.	+	.	gene_id "ENSG00000223972"; transcript_id "ENST00000450305"; gene_name "DDX11L1"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "DDX11L1-001"; transcript_source "havana";
1	unprocessed_pseudogene	transcript	14363	29370	.	-	.	gene_id "ENSG00000227232"; transcript_id "ENST00000438504"; gene_name "WASH7P"; gene_source "ensembl_havana"; gene_biotype "pseudogene"; transcript_name "WASH7P-202"; transcript_source "ensembl";


Ensembl hg38: 

0			1	       2		3	    4     5	    6 7		
1       havana  transcript    11869   14409   .     + .       

  8_0		 8_1				8_2		8_3		8_4		            8_5		         8_6		 8_7	 8_8	  8_9         8_10	     8_11	     8_12                8_13							8_14		   8_15				8_16			8_17	   8_18					8_19				
gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; 

tag "basic"; transcript_support_level "1";

Arg 2 - Op file name
Arg 3 - Whether it is a reference ensembl annotation or not  -- can be either "ens", "ens38" or "gen" to indicate ensembl or gencode resp. or any other word to indicate not ensembl/gencode
"""

tid2gid = []
mgid2hgnc = {}
mgid2egid = {}
mgid2coords = {}

with open(sys.argv[1], "r") as inph:
	for line in inph:
		if "#" in line:
			continue
		rec = line.strip().split("\t")
		if rec[2] != "transcript":
			continue

		#print(rec)
		rec8 = rec[8].replace(";", "").replace('"', '').split(" ")
		
		mg = rec8[1]
		if mg in mgid2coords:
			a, b = mgid2coords[mg]
			mgid2coords[mg] = (min(int(rec[3]), b), max(int(rec[4]), b))
		else:
			mgid2coords[mg] = (int(rec[3]), int(rec[4]))

		if sys.argv[3] not in ["ens", "ens38", "gen"]:
			tid2gid.append([rec[0].replace("chr", ""), rec[6], rec8[3], rec8[1]])

			if len(rec8) > 4:
				if mg not in mgid2egid:
					mgid2egid[mg] = []
		
				if mg not in mgid2hgnc:
					mgid2hgnc[mg] = []
		
				mgid2egid[mg].append(rec8[7])
				mgid2hgnc[mg].append(rec8[5])
			continue
		
		if mg not in mgid2egid:
			mgid2egid[mg] = []
		
		if mg not in mgid2hgnc:
			mgid2hgnc[mg] = []
		
		#print(rec8)
		#sys.exit()
		if sys.argv[3] == "gen":
			tid2gid.append([rec[0].replace("chr", ""), rec[6], rec8[3], rec8[1], rec8[9], rec8[5]])
			mgid2egid[mg].append(rec8[1])
			mgid2hgnc[mg].append(rec8[7])
			continue

		if sys.argv[3] == "ens":
			tid2gid.append([rec[0].replace("chr", ""), rec[6], rec8[3], rec8[1], rec[1], rec8[9]])
			mgid2egid[mg].append(rec8[1])
			mgid2hgnc[mg].append(rec8[5])
			continue
		
		if sys.argv[3] == "ens38":
			tid2gid.append([rec[0].replace("chr", ""), rec[6], rec8[5], rec8[1], rec8[19], rec8[13]])
			mgid2egid[mg].append(rec8[1])
			mgid2hgnc[mg].append(rec8[9])
			continue

with open(sys.argv[2], "w") as outh:
	if sys.argv[3] in ["ens", "ens38", "gen"]:
		outh.write("TID\tGID\tCHROM\tSTRAND\tSTART\tEND\tEGID\tHGNC\tTXBIO\tGNBIOT\n")
		for tg in tid2gid:
			#print(tg)
			#sys.exit()
			ch = tg[0]
			strand = tg[1]
			t = tg[2]
			mg = tg[3]
			eg = ",".join(list(set(mgid2egid[mg])))
			hg = ",".join(list(set(mgid2hgnc[mg])))
			txb = tg[-2]
			gnb = tg[-1]

			s, e = mgid2coords[mg][0], mgid2coords[mg][1] 
			
			op = t+"\t"+mg+"\t"+ch+"\t"+strand+"\t"+str(s)+"\t"+str(e)+"\t"+eg+"\t"+hg+"\t"+txb+"\t"+gnb
			outh.write(op+"\n")

	else:
		outh.write("TID\tGID\tCHROM\tSTRAND\tSTART\tEND\tEGID\tHGNC\n")
		for tg in tid2gid:
			ch = tg[0]
			strand = tg[1]
			t = tg[2]
			mg = tg[3]
			eg = tg[3]
			hg = tg[3]
			s, e = mgid2coords[mg][0], mgid2coords[mg][1] 
			
			if mg in mgid2egid:
				eg = ",".join(list(set(mgid2egid[mg])))
			
			if mg in mgid2hgnc:
				hg = ",".join(list(set(mgid2hgnc[mg])))
				
			op = t+"\t"+mg+"\t"+ch+"\t"+strand+"\t"+str(s)+"\t"+str(e)+"\t"+eg+"\t"+hg
			outh.write(op+"\n")
			
