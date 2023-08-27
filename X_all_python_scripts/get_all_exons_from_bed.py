import sys

"""
1	14362	29370	ENST00000423562	0	-	14362	29370	0	10	467,69,152,159,198,136,137,147,154,50,	0,607,1433,2244,2495,2870,3243,3552,10375,14958,
1	14466	29331	MSTRG.2.6	0	-	14466	29331	0	11	363,69,152,159,202,136,137,147,112,154,11,	0,503,1329,2140,2387,2766,3139,3448,3801,10271,14854,
1	96644	105848	MSTRG.5.1	0	+	96644	105848	0	2	166,1213,	0,7991,
1	134900	139379	ENST00000423372	0	-	134900	139379	0	2	902,1759,	0,2720,
1	146385	173862	ENST00000466557	0	-	146385	173862	0	8	124,65,529,59,66,216,132,110,	0,9381,17877,19498,21714,22663,26171,27367,
"""

exonsD = {}

with open(sys.argv[2], "w") as oph:
	with open(sys.argv[1], "r") as inph:
		for line in inph:
			rec = line.strip().split("\t")
			s, e = int(rec[1]), int(rec[2])
			exons = [(s + int(es), s + int(es) + int(el)) for (es, el) in zip(rec[11].split(",")[:-1], rec[10].split(",")[:-1])]

			for ex in exons:
				idx = "__".join([str(e) for e in [rec[0], ex[0], ex[1], rec[5]]])

				if idx not in exonsD:
					exonsD[idx] = [rec[0], ex[0], ex[1], rec[5]]
	exid = 1
	for idx in exonsD:
		orec = [exonsD[idx][0], exonsD[idx][1], exonsD[idx][2], "exon."+str(exid), 0, exonsD[idx][-1]] 
		oph.write("\t".join([str(e) for e in orec]) + "\n")

		exid += 1

