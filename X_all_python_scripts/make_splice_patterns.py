import sys
import re
from collections import Counter

"""
 0   1		2		3			4		 5		 6	        7	 8		9         10		   11			12		  13		14	15
EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
D3.event2	1	+	764484	787307	D3.event2.1,D3.event2.2,D3.event2.3	,783034-783186^,783034-783186^784864-784982^	MSTRG.93	ENSG00000228794,ENSG00000234711,ENSG00000269308	LINC01128,TUBB8P11,AL645608.2	764484	787307	D3.event2.1,D3.event2.2,D3.event2.3	,783034-783186^,783034-783186^784864-784982^	3	3
D2.event3	1	+	787490	788771	D2.event3.1,D2.event3.2	,788051-788146^	MSTRG.93	ENSG00000228794,ENSG00000234711,ENSG00000269308	LINC01128,TUBB8P11,AL645608.2	787490	788771	D2.event3.1,D2.event3.2	,788051-788146^	2	2
"""

other = {"^":"-", "-":"^"}

def get_splice_pattern(splc, strand):
	
	splcs = splc.split(",")  ## [6874384-6874198^6873621-6873536^6856695-6856645^6855309, 6855285-]
	final = []
	all_sites = []

	##first add any missing symbols
	for i in range(len(splcs)):   ## len = 2
		sc = splcs[i] 
		if sc != "":
			if sc[-1] not in ["^", "-"]:
				penultimate = re.sub(r'\d', '', sc)[-1] ## -^-^-^ >> penult = ^ >> other[^] = -
				sc = sc + other[penultimate] # 6874384-6874198^6873621-6873536^6856695-6856645^6855309 + "-" >> 6874384-6874198^6873621-6873536^6856695-6856645^6855309-
		
			splcs[i] = sc ## 6874384-6874198^6873621-6873536^6856695-6856645^6855309-
			all_sites.extend([int(e) for e in re.sub(r'[\^\-]', ',', sc)[:-1].split(",")]) ## 6874384,6874198,6873621,6873536,6856695,6856645,6855309, >> [6874384,6874198,6873621,6873536,6856695,6856645,6855309]

	all_sites = sorted(list(set(all_sites)))
	if strand == "-":
		all_sites = all_sites[::-1]
	for sc in splcs:
		f = ""
		if sc != "":
			spl_sites = [all_sites.index(int(e))+1 for e in re.sub(r'[\^\-]', ',', sc)[:-1].split(",")] ## 
			for (g,s) in zip(spl_sites, re.sub(r'\d', '', sc)):
				f += str(g)+s
			
		final.append(f)
	
	return ",".join(sorted(final))

"""
example
splc = "6874384-6874198^6873621-6873536^6856695-6856645^6855309,6855285-,"
print(get_splice_pattern(splc, "-"))
"""

ev = []

with open(sys.argv[2] + ".symb.splc.out", "w") as outh:
	with open(sys.argv[1], "r") as inph:
		for line in inph:
			if "EID" in line:
				outh.write(line.strip()+"\tSYMB.SPLC\n")
				continue
			
			rec = line.strip().split("\t")
			strand, splc = rec[2], rec[13]

			symbsplc = get_splice_pattern(splc, strand)
			
			outh.write("\t".join(rec) + "\t" + symbsplc + "\n")

			ev.append(symbsplc)


ev = Counter(ev)

with open(sys.argv[2] + ".symb.freq.out", "w") as outh:
	outh.write("SYMB.SPLC\tFREQUENCY\tDIM\n")
	for sym in ev:
		dim = len(sym.split(","))
		
		outh.write(sym + "\t" + str(ev[sym]) + "\t" + str(dim) + "\n")

