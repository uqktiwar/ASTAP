import sys
from operator import itemgetter
import re

"""
 0	 1		2		3			4		 5			6		 7	 8		  9			10			11			12			13		14	15
EID	CHR	STRAND	OFLANK.5P	OFLANK.3P	ISOS	SPLCHAIN	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD
D2.event1	1	+	783186	787307	D2.event1.1,D2.event1.2	,784864-784982^	MSTRG.93	ENSG00000228794,ENSG00000234711,ENSG00000269308	LINC01128,TUBB8P11,AL645608.2	783186	787307	D2.event1.1,D2.event1.2	,784864-784982^	2	2
D3.event2	1	+	764484	787307	D3.event2.1,D3.event2.2,D3.event2.3	,783034-783186^,783034-783186^784864-784982^	MSTRG.93	ENSG00000228794,ENSG00000234711,ENSG00000269308	LINC01128,TUBB8P11,AL645608.2	764484	787307	D3.event2.1,D3.event2.2,D3.event2.3	,783034-783186^,783034-783186^784864-784982^	3	3

D6.event21646	X	-	71694579	71791907	
D6.event21646.1,D6.event21646.2,D6.event21646.3,D6.event21646.4,D6.event21646.5,D6.event21646.6	
71788734-71788604^71787880-71787739^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71700300-71700300^,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710775^71709113-71709112^71708885-71708783^71700300-71700300^,71788734-71788604^71787880-71787739^71715118-71715006^71708891-71708783^71700300-71700300^,71715118-71715006^71708891-71708783^	
MSTRG.18776	
ENSG00000147099	
HDAC8	71700300	71791907	
D6.event21646.1,D6.event21646.2,D6.event21646.3,D6.event21646.4,D6.event21646.5,D6.event21646.6	
71788734-71788604^71787880-71787739^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710779^71708891-71708783^71700300-71700300^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71710856-71710775^71709113-71709112^71708885-71708783^71700300-71700300^71694579f,71788734-71788604^71787880-71787739^71715118-71715006^71708891-71708783^71700300-71700300^71694579f,71715118-71715006^71708891-71708783^71694579f	6	6

"""


def get_chrwise_dimwise_events(e2gpath):

	events = {}
	headr = None
	maxdim = 0
	with open(e2gpath, 'r') as inph:
		for line in inph:

			rec = line.strip().split("\t")

			if rec[0] == "EID":
				headr = line # store header for the final output file
				continue
			
			if rec[1] not in events:
				events[rec[1]] = {"+":{}, "-":{}}

			rd = int(rec[15])

			maxdim = max(maxdim, rd)
			
			if rd not in events[rec[1]][rec[2]]:
				events[rec[1]][rec[2]][rd] = []
			
			events[rec[1]][rec[2]][rd].append([rec[0], int(rec[10]), int(rec[11]), rec[13]])
	
	return (headr, maxdim, events)


def is_substring_of_splchain(iso, splchain):
	
	for splc in splchain:
		if iso in splc:
			###print(iso, "found")
			return True
	
	return False

def is_contained_in_splchain(f02s, f12s, splchain2, f01s, f11s, splchain1):
		
	splchain1 = [f01s+sc+f11s for sc in splchain1.split(",")]
	##print(splchain1)

	for iso in splchain2.split(","):
		iso = f02s+iso+f12s
		###print(iso, "checking")
		if (iso not in splchain1) and (is_substring_of_splchain(iso, splchain1) == False):
			return False
	
	return True


def append_symb_to_flanks(strand, event):
	fsyms = {"^":"-", "-":"^"}
	##print(event)
	eid, f0, f1, splchain = event
	
	splchain = re.sub(r'[\d\,]', '', splchain).replace("[", "-").replace("]", "^")
	
	if strand == "+":
		return (str(f0)+fsyms[splchain[0]], str(f1)+fsyms[splchain[-1]])
	
	return (str(f1)+fsyms[splchain[0]], str(f0)+fsyms[splchain[-1]])
		


e2gpath = sys.argv[1]

headr, maxdim, events = get_chrwise_dimwise_events(e2gpath)

rmkeys = []

for chrom in events:
	for strand in events[chrom]:
		for d1 in range(maxdim, 1, -1):
			if d1 in events[chrom][strand]:
				events[chrom][strand][d1] = sorted(events[chrom][strand][d1], key = itemgetter(1, 2), reverse = False)


for chrom in events:
	for strand in events[chrom]:
		for d1 in range(maxdim, 1, -1):
			if d1 in events[chrom][strand]:
				for event1 in events[chrom][strand][d1]:
					eid1, f01, f11, sp1 = event1
					f01s, f11s = append_symb_to_flanks(strand, event1) ## IMPORTANT !!!! here f0s and f1s are returned in the order they need to be appended to the splchain .. i.e. for events on neg strand numeric part of (f1s) < numeric part of (f0s)

					if eid1 not in rmkeys:
						for d2 in range(d1, 1, -1):
							if d2 in events[chrom][strand]:
								for event2 in events[chrom][strand][d2]:
									eid2, f02, f12, sp2 = event2
									f02s, f12s = append_symb_to_flanks(strand, event2) ## IMPORTANT !!!! here f0s and f1s are returned in the order they need to be appended to the splchain .. i.e. for events on neg strand numeric part of (f1s) < numeric part of (f0s)
									
									##print(event1, event2)

									if eid1 != eid2:
										if f02 >= f11:
											break

										if f02 >= f01 and f12 <= f11:
											###print(is_contained_in_splchain(f02s, f12s, sp2, f01s, f11s, sp1))
											if is_contained_in_splchain(f02s, f12s, sp2, f01s, f11s, sp1):
												rmkeys.append(eid2)


op1 = open(sys.argv[2]+".final.e2g.out", "w")
op1.write(headr)

op2 = open(sys.argv[2]+".events2requant.modFlanks.out", "w")
op2.write(headr)

with open(e2gpath, "r") as inph:
	for line in inph:

		rec = line.strip().split("\t")

		if rec[0] == "EID":
			continue
			
		if rec[0] not in rmkeys:
			op1.write(line)
		
		if (rec[10], rec[11]) != (rec[3], rec[4]):
			op2.write(line)


op1.close()
op2.close()

"""								
event1 = ["D4.event7028", 438195, 442653, "441104-440985^,441104-440960^,441096-440985^"]
eid1, f01, f11, sp1 = ["D4.event7028", 438195, 442653, "441104-440985^,441104-440960^,441096-440985^"]
f01s, f11s = append_symb_to_flanks("-", event1) ## IMPORTANT !!!! here f0s and f1s are returned in the order they need to be appended to the splchain .. i.e. for events on neg strand numeric part of (f1s) < numeric part of (f0s)

#print(f01s, f11s)

event2 = ["D3.event7027", 438195, 442653,"441104-440985^,441104-440960^,441096-440985^"]
eid2, f02, f12, sp2 = event2
f02s, f12s = append_symb_to_flanks("-", event2) ## IMPORTANT !!!! here f0s and f1s are returned in the order they need to be appended to the splchain .. i.e. for events on neg strand numeric part of (f1s) < numeric part of (f0s)

#print(f02s, f12s)

if f02 >= f01 and f12 <= f11:
	#print(is_contained_in_splchain(f02s, f12s, sp2, f01s, f11s, sp1))


D3.event7027	12	-	438195	442653	D3.event7027.1,D3.event7027.2,D3.event7027.3	441104-440985^,441104-440960^,441096-440985^	MSTRG.33613	ENSG00000073614	KDM5A	438195	442653	D3.event7027.1,D3.event7027.2,D3.event7027.3	441104-440985^,441104-440960^,441096-440985^	3	3
	D4.event7028.1,D4.event7028.2,D4.event7028.3,D4.event7028.4	442815-442653^441104-440985^,442815-442653^441104-440960^,442815-442653^441096-440985^,441096-440985^	MSTRG.33613	ENSG00000073614	KDM5A	438195	442653	D4.event7028.1,D4.event7028.2,D4.event7028.3	441104-440985^,441104-440960^,441096-440985^	4	3
D2.event7029	12	-	417171	418969	D2.event7029.1,D2.event7029.2	,418728-418686^	MSTRG.33613	ENSG00000073614	KDM5A	417171	418969	D2.event7029.1,D2.event7029.2	,418728-418686^	2	2
D3.event7030	12	-	416620	420051	D3.event7030.1,D3.event7030.2,D3.event7030.3	419130-418969^418728-418686^417171-,419130-418969^417171-,417021-	MSTRG.33613	ENSG00000073614	KDM5A	417171	418969	D3.event7030.1,D3.event7030.2	418728-418686^,	3	2
D2.event7031	12	-	404739	406207	D2.event7031.1,D2.event7031.2	404974-,404959-	MSTRG.33613	ENSG00000073614	KDM5A	404739	406207	D2.event7031.1,D2.event7031.2	404974-,404959-	2	2
D2.event7032	12	-	394828	401925	D2.event7032.1,D2.event7032.2	,395373-395292^	MSTRG.33613	ENSG00000073614	KDM5A	394828	401925	D2.event7032.1,D2.event7032.2	,395373-395292^	2	2
D2.event7033	12	-	1036429	1039052	D2.event7033.1,D2.event7033.2	1038985^,1038930^	MSTRG.33659	ENSG00000250132,ENSG00000002016	RAD52,RP11-359B12.2	1036429	1039052	D2.event7033.1,D2.event7033.2	1038985^,1038930^	2	2
D2.event7034	12	-	1023287	1025510	D2.event7034.1,D2.event7034.2	,1023698-1023597^	MSTRG.33659	ENSG00000250132,ENSG00000002016	RAD52,RP11-359B12.2	1023287	1025510	D2.event7034.1,D2.event7034.2	,1023698-1023597^	2	2
D3.event7037	12	-	2787063	2799306	D3.event7037.1,D3.event7037.2,D3.event7037.3	,2791778-2791652^2791014-2790820^2787304-2787185^,2787304-2787185^	MSTRG.33822	ENSG00000256271,ENSG00000246627	CACNA1C-AS1,CACNA1C-AS2	2787063	2799306	D3.event7037.1,D3.event7037.2,D3.event7037.3	,2791778-2791652^2791014-2790820^2787304-2787185^,2787304-2787185^	3	3
D4.event7038	12	-	2787063	2799376	D4.event7038.1,D4.event7038.2,D4.event7038.3,D4.event7038.4	2799306^,2799306^2791778-2791652^2791014-2790820^2787304-2787185^,2799306^2787304-2787185^,2799243^	MSTRG.33822	ENSG00000256271,ENSG00000246627	CACNA1C-AS1,CACNA1C-AS2	2787063	2799306	D4.event7038.1,D4.event7038.2,D4.event7038.3	,2791778-2791652^2791014-2790820^2787304-2787185^,2787304-2787185^	4	3
"""


