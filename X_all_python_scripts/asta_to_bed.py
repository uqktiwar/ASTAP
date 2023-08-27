import sys
import re
import  math 

"""
Arg1 - asta gtf file
0		1			2		  3		   4	5	6	7		8_0								8_1					   8_2			8_3				 8_4		   8_5			 8_6	   8_7	       8_8				8_9			   8_10			8_11			   8_12	   8_13	  8_14	 8_15	8_16	8_17	
1	Undefined	as_event	7837220	7838229	.	+	.	transcript_id "MSTRG.462.3,ENST00000054666/MSTRG.462.1"; gene_id "1:7831329-7841492W"; flanks "7837220-,7838229^"; structure "0,1^2-"; splice_chain ",7837378^7838178-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 
1	Undefined	as_event	8044954	8045339	.	+	.	transcript_id "ENST00000338639/MSTRG.517.17,MSTRG.517.11"; gene_id "1:8014351-8429601W"; flanks "8044954-,8045339]"; structure "0,1^2-"; splice_chain ",8044976^8045151-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 
1	Undefined	as_event	9781185	9781645	.	+	.	transcript_id "MSTRG.712.7,MSTRG.712.18/MSTRG.712.3"; gene_id "1:9689987-9788977W"; flanks "9781185-,9781645^"; structure "0,1^2-"; splice_chain ",9781306^9781502-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 
1	Undefined	as_event	10228194	10231387	.	+	.	transcript_id "MSTRG.744.19,ENST00000377157/MSTRG.744.18"; gene_id "1:10093016-10241298W"; flanks "10228194-,10231387^"; structure "0,1^2-"; splice_chain ",10228328^10231196-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 
1	Undefined	as_event	10479474	10480201	.	+	.	transcript_id "MSTRG.775.1,ENST00000270776"; gene_id "1:10448367-10480568W"; flanks "10479474-,10480201]"; structure "0,1^2-"; splice_chain ",10479596^10479687-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 
1	Undefined	as_event	11080486	11083386	.	+	.	transcript_id "ENST00000473869/MSTRG.584.9,MSTRG.584.14"; gene_id "1:11072699-11090585W"; flanks "11080486-,11083386^"; structure "1^4-5^6-,2^3-"; splice_chain "11080656^11082181-11082307^11083250-,11080684^11081972-"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 

1	Undefined	as_event	796653	817432	.	+	.	transcript_id "MSTRG.93.20,MSTRG.93.16"; gene_id "1:564442-832582W"; flanks "796653^,817432-"; structure "0,1-2^"; splice_chain ",800640-800879^"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 

Arg2 - op filename
"""


def get_isowise_exons(rec8, strand):
	
	f0, f1 = rec8[5].split(",")  ## 804055-, 800637^
	f0, f1 = f0.replace("[", "-").replace("]", "^"), f1.replace("[", "-").replace("]", "^") 

	f00, f11 = "", ""
	print(rec8)
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

	splc = rec8[9].split(",") ## ["", "802576^802469-"]

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
		eid = 0
		for line in inph:
			eid += 1
			print(eid)
			rec = line.strip().split("\t")

			rec8 = rec[8].replace(";", "").replace('"', '').split(" ")
			dim = int(math.sqrt(int(rec8[15])))

			support = [len(txs.split("/")) for txs in rec8[1].split(",")]   ## to be used as a support score for the different isoforms

			exons = get_isowise_exons(rec8, rec[6])

			s, e = min([ex[0] for ex in exons[0]]) - 1, max([ex[1] for ex in exons[0]])
			
			## 11	75112683	75115251	event3.4	0	+	75112683	75115251	0	3	94,95,188,	0,712,2380,
			tid = 0
			for i in range(len(exons)):
				txex = exons[i]
				estarts = [str(ex[0]-1-s) for ex in txex]
				elens = [str(1+ex[1]-ex[0]) for ex in txex]

				tid += 1
				event = "D" + str(dim) + ".event" + str(eid) + "." + str(tid) 
				orec = [rec[0], str(s), str(e), event, str(support[i]), rec[6], str(s), str(e), str(support[i]), str(len(txex)), ",".join(elens)+",", ",".join(estarts)+","]
				
				outh.write("\t".join(orec)+"\n")




			
			
			

			












