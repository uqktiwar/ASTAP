import sys
import os

"""
0	  1			2		3		4	5	  6			7	8	9	  10	  11
1	800639	802849	D2.event1.1	1	+	800639	802849	1	2	240,50,	0,2160,
1	800639	802849	D2.event1.2	1	+	800639	802849	1	2	590,50,	0,2160,
1	796603	817481	D2.event2.1	1	+	796603	817481	1	2	50,50,	0,20828,
1	796603	817481	D2.event2.2	1	+	796603	817481	1	3	50,240,50,	0,4036,20828,
"""


inrs = {}
with open(os.path.join(sys.argv[1], sys.argv[2]), "r") as inph:
	for line in inph:

		rec = line.strip().split("\t")
		s, e = int(rec[1]), int(rec[2])

		estarts = [int(es) for es in rec[11].split(",")[:-1]]
		elens = [int(el) for el in rec[10].split(",")[:-1]]

		exons = [(s+es, s+es+el) for (es, el) in zip(estarts, elens)]

		introns = [(exons[i][1], exons[i+1][0]) for i in range(len(exons)-1)]
		
		for inr in introns:
			inrdx = "__".join([rec[0], str(inr[0]), str(inr[1]), rec[5]])

			if inrdx not in inrs:
				inrs[inrdx] = [rec[0], str(inr[0]), str(inr[1]), rec[4], rec[5]]


inrdx = 0
with open(os.path.join(sys.argv[1], sys.argv[3]+".IR.junctions.bed"), "w") as outh:
	for inr in inrs:
		inrdx += 1

		orec = inrs[inr][0:3] + ["intron."+str(inrdx)] + inrs[inr][3:] 
		outh.write("\t".join(orec)+"\n")





