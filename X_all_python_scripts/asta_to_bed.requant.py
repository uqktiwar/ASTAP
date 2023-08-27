import sys
import re
import  math 

"""
Arg1 - asta gtf file
 0	 1	   2		3			4		  5		   6		 7	  8		 9			10			11			12			13		14	15
EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
D3.event167	1	+	31769922	31791956	D3.event167.1,D3.event167.2,D3.event167.3	,31782891-31783011^,31782891-31783011^31791014-31791083^	MSTRG.649	ENSG00000121766	ZCCHC17	31783011	31791956	D3.event167.2,D3.event167.3	,31791014-31791083^	3	2
D3.event178	1	+	32688231	32694100	D3.event178.1,D3.event178.2,D3.event178.3	32689635-32689722^32690011-32690076^,32689635-32689722^32690011-32690076^32691772-32691921^32692004-32692131^,32691772-32691921^32692004-32692131^	MSTRG.676	ENSG00000084623	EIF3I	32690076	32694100	D3.event178.1,D3.event178.2	,32691772-32691921^32692004-32692131^	3	2
D3.event220	1	+	39457272	39469019	D3.event220.1,D3.event220.2,D3.event220.3	39463843-39463983^,39463843-39463983^39466644-39466778^,39466644-39466778^	MSTRG.796	ENSG00000174574	AKIRIN1	39463983	39469019	D3.event220.1,D3.event220.2	,39466644-39466778^	3	2

Arg2 - op filename
"""

def add_flank_symbols(f0, f1, splc, strand):

	flank2sym = {"^":"-", "-":"^"}

	splsym = re.sub(r'[\d\,]', '', splc)
	splsym = splsym.replace("[", "-").replace("]", "^")

	f0s, f1s = flank2sym[splsym[0]], flank2sym[splsym[-1]]
	
	if strand == "+":
		return (f0+f0s, f1+f1s)
	
	if strand == "-":
		return (f1+f0s, f0+f1s)


def get_isowise_exons(rec, strand):
	
	f0, f1 = rec[10], rec[11]
	f0, f1 = add_flank_symbols(f0, f1, splc = rec[13], strand = rec[2])

	f00, f11 = "", ""
	##print(rec)
	if f0[-1] == "^":
		if strand == "+":
			f00 = str(int(f0[:-1]) - 49) + "-"
		if strand == "-":
			f00 = str(int(f0[:-1]) + 49) + "-"
	
	if f1[-1] == "-":
		if strand == "+":
			f11 = str(int(f1[:-1]) + 49) + "^"
		if strand == "-":
			f11 = str(int(f1[:-1]) - 49) + "^"
	
	f0, f1 = f00+f0, f1+f11

	splc = rec[13].split(",") ## ["", "802576^802469-"]

	pex = re.compile(r'\d+\-\d+\^')
	exons = []
	
	for sc in splc:
		scf = f0+sc+f1

		print(scf) ## [804055-800637^, 804055-802576^802469-800637^]
		
		txex = [m.group() for m in pex.finditer(scf)] ## [804055-802576^, 802469-800637^] >> [804055,802576,, 802469,800637,]
		txex = [[int(exi) for exi in re.sub(r'[\^\-]', ',', ex).split(",")[:-1]] for ex in txex] ##[[804055,802576], [802469,800637]]
		txex = sorted([(min(ex), max(ex)) for ex in txex]) ##[(802576, 804055), (800637, 802469)]

		exons.append(txex)
	
	return exons

with open(sys.argv[2]+".bed", "w") as outh:
	with open(sys.argv[1], "r") as inph:
		#eid = 0
		for line in inph:
			rec = line.strip().split("\t")
			if rec[0] == "EID":
				continue

			eid = rec[0]
			print(eid)
			dim = int(rec[15])

			support = 1 ## Can't do this --> [len(txs.split("/")) for txs in rec8[1].split(",")] as this info is not available in the new splice/e2g format file 

			exons = get_isowise_exons(rec, rec[2])

			s, e = min([ex[0] for ex in exons[0]]) - 1, max([ex[1] for ex in exons[0]])
			
			## 11	75112683	75115251	event3.4	0	+	75112683	75115251	0	3	94,95,188,	0,712,2380,
			for i in range(len(exons)):
				txex = exons[i]
				estarts = [str(ex[0]-1-s) for ex in txex]
				elens = [str(1+ex[1]-ex[0]) for ex in txex]

				tid = rec[12].split(",")[i]
				event = tid 
				orec = [rec[1], str(s), str(e), event, str(support), rec[2], str(s), str(e), str(support), str(len(txex)), ",".join(elens)+",", ",".join(estarts)+","]
				
				outh.write("\t".join(orec)+"\n")




			
			
			

			












