import sys

"""
Sample input
chr1	Undefined	as_event	1424584	1426051	.	+	.	transcript_id "ENST00000474481,ENST00000485748,ENST00000308647/MSTRG.70.12/MSTRG.70.14"; gene_id "chr1:1390840-1470065W"; flanks "1424584-,1426051^"; structure "0,1^2-3^4-5^6-,1^4-5^6-"; splice_chain ",1424654^1425072-1425191^1425637-1425804^1425943-,1424654^1425637-1425804^1425943-"; sources "Undefined,Undefined,Undefined"; NMorMRNA "null"; degree "9"; dimension "3_3"; 
 0			1			2		   3	   4	5	6	7	    8_0													8_1									8_2			  8_3				8_4			  8_5			8_6					8_7					8_8											8_9

sys.argv[1] -- Event file in ASTA GTF format
sys.argv[2] -- Tx2gene table for the transcripts used to mine the AS events Format: "TID\tGID\tCHROM\tSTRAND\tSTART\tEND\tEGID\tHGNC\n"
sys.argv[3] -- Name of the output event 2 gene mapping file
"""

tx2g = {}

with open(sys.argv[2], "r") as tx2gf:
	for line in tx2gf:
		## "TID\tGID\tCHROM\tSTRAND\tSTART\tEND\tEGID\tHGNC\n"
		rec = line.strip().split("\t")

		tx2g[rec[0]] = [rec[1], rec[6], rec[7]]
						#GID	#EGID	#HGNC


eid = 1
with open(sys.argv[3], "w") as outh:
	outh.write("EID\tCHR\tSTRAND\tFLANK.5P\tFLANK.3P\tSPLCHAIN\tGID\tEGID\tHGNC\n")
	with open(sys.argv[1], "r") as inph:
		eid = 1
		for line in inph:
			rec = line.strip().split("\t")
			rec8 = rec[8].replace(";", "").replace('"', '').split(" ")
			
			print(rec8)

			t = rec8[1].split(",")[-1].split("/")[0]
			txRecs = tx2g[t]

			op = "\t".join(["event"+str(eid), rec[0], rec[6], rec[3], rec[4], rec8[9]] + txRecs) + "\n"
			outh.write(op)

			eid += 1

