import sys
import re

"""
 0	 1	   2		3			4		 5			6		 7	  8		 9			10			11			12			13		14	15
EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
D2.event1	1	+	321032	322038	D2.event1.1,D2.event1.2	321264^,321290^	MSTRG.17	ENSG00000237094,ENSG00000250575,ENSG00000224813	RP4-669L17.10,RP4-669L17.4,RP4-669L17.8	321032	322038	D2.event1.1,D2.event1.2	321264^,321290^	2	2

"""


eids = []

with open(sys.argv[2]+".microex.txt", "w") as oph:
	oph.write("EID\tMICRO\n")
	with open(sys.argv[1], "r") as inph:
		for line in inph:

			if "STRAND" in line:
				continue
			
			rec = line.strip().split("\t")

			splc = rec[13]

			micro = 0

			for sc in splc.split(","):
				for m in re.finditer(r'\d+\-\d+\^', sc):
					exm = [int(e) for e in m.group()[:-1].split("-")]
					exm = sorted(exm)

					elen = exm[1] - exm[0] + 1
					
					if elen <= 50:
						micro = 1
						break	
				if micro:
					break
			

			oph.write(rec[0]+"\t"+str(micro)+"\n")
		
				
