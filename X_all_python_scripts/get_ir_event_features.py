import sys
import os

"""
0		1		2		3			4	5		6		7	8	9		10		 11			12
1	1148371	1149165	D2.event739.2	0	-	1148371	1149165	0	2	102,123,	0,671,	D2.event739
1	1148371	1149165	D2.event739.1	0	-	1148371	1149165	0	1	794,	0,	D2.event739
1	1682670	1684499	D2.event760.1	0	-	1682670	1684499	0	1	1829,	0,	D2.event760
1	1682670	1684499	D2.event760.2	0	-	1682670	1684499	0	2	1166,619,	0,1210,	D2.event760
1	3541565	3541837	D2.event10.1	0	+	3541565	3541837	0	1	272,	0,	D2.event10
1	3541565	3541837	D2.event10.2	0	+	3541565	3541837	0	2	73,82,	0,190,	D2.event10
"""

def map_iso_to_event(bedpath):

	iso2ev = {}

	with open(bedpath, "r") as bedh:
		for line in bedh:

			rec = line.strip().split("\t")
			if len(rec) == 12:
				rec.append(".".join(rec[3].split(".")[:-1]))

			if rec[-1] not in iso2ev:
				iso2ev[rec[-1]] = {"chrom":rec[0], "estart":rec[1], "eend":rec[2], "strand":rec[5], "exons":[], "introns":[], "iso":{}}
			
			s, e = int(rec[1]), int(rec[2])
			exons = [(s + int(es), s + int(es) + int(el)) for (es, el) in zip(rec[11].split(",")[:-1], rec[10].split(",")[:-1])]
			introns = [(exons[i][1], exons[i+1][0]) for i in range(len(exons) - 1)]

			iso2ev[rec[-1]]["exons"].extend(exons)
			iso2ev[rec[-1]]["introns"].extend(introns)
			
			iso2ev[rec[-1]]["iso"][rec[3]] = {"exons":exons, "introns":introns, "ret":[]}
	
	for eid in iso2ev:
		iso2ev[rec[-1]]["exons"] = list(set(iso2ev[rec[-1]]["exons"]))
		iso2ev[rec[-1]]["introns"] = list(set(iso2ev[rec[-1]]["introns"]))

	
	return iso2ev


def find_retained_introns(iso2ev):
	
	for eid in iso2ev:
		for inr in iso2ev[eid]["introns"]:
			for iso in iso2ev[eid]["iso"]:
				for ex in iso2ev[eid]["iso"][iso]["exons"]:
					if ex[0] < inr[0] and inr[1] < ex[1]:
						iso2ev[eid]["iso"][iso]["ret"].append(inr)
	

	return iso2ev

iso2ev = map_iso_to_event(os.path.join(sys.argv[1], sys.argv[2]))
iso2ev = find_retained_introns(iso2ev)

with open(os.path.join(sys.argv[1], "all_retained_introns.bed"), "w") as irh:
	for eid in iso2ev:
		for iso in iso2ev[eid]["iso"]:
			for irt in iso2ev[eid]["iso"][iso]["ret"]:
				orec = [iso2ev[eid]["chrom"], irt[0], irt[1], iso, 0, iso2ev[eid]["strand"]]

				irh.write("\t".join([str(e) for e in orec]) + "\n")


with open(os.path.join(sys.argv[1], "all_ir_events_features.out"), "w") as oph:
	oph.write("\t".join(["gene_id", "feature_id", "chrom", "strand", "ev.start", "ev.end", "exons", "junctions", "retained_introns"]) + "\n")
	for eid in iso2ev:
		for iso in iso2ev[eid]["iso"]:
			exons = ",".join([str(e) for e in iso2ev[eid]["iso"][iso]["exons"]])
			introns = ",".join([str(e) for e in iso2ev[eid]["iso"][iso]["introns"]])
			ret = ",".join([str(e) for e in iso2ev[eid]["iso"][iso]["ret"]])
			orec = [eid, iso, iso2ev[eid]["chrom"], iso2ev[eid]["strand"], iso2ev[eid]["estart"], iso2ev[eid]["eend"], exons, introns, ret]

			oph.write("\t".join([str(e) for e in orec]) + "\n")










			