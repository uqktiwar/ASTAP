import sys
from collections import Counter
import re
"""
0	  1		  2				3		4	5	   6	   7	8	9	   10				11
1	11868	14409	ENST00000456328	0	+	11868	14409	0	3	359,109,1189,	0,744,1352,
1	11871	14412	ENST00000515242	0	+	11871	14412	0	3	356,109,1188,	0,741,1353,
"""

def map_exons_to_tx(txbed):

	exons_to_tx = {}

	with open(txbed, "r") as bedh:
		for line in bedh:

			rec = line.strip().split("\t")
			
			e, s = int(rec[1]), int(rec[2])

			exons = [(e + int(es) + 1, e + int(es) + int(el)) for (es, el) in zip(rec[11].split(",")[:-1], rec[10].split(",")[:-1])]

			exons_to_tx[rec[3]] = {"chrom":rec[0], "strand":rec[5], "exons":exons}
	
	return exons_to_tx


def find_best_flanking_exon(f0, f1, supp_txs, exons_to_tx):
	"""
	f0, f1 - integral flanks low (f0) and high (f1) 
	supp_txs - list of all transcript names supporting the event [tx1, tx2, ..]
	exons_to_tx - exons to tx mapping generated using map_exons_to_tx function
	"""
	f0, f1 = int(f0), int(f1)
	f0exons, f1exons = [], []
	for tx in supp_txs:
		for exon in exons_to_tx[tx]["exons"]:

			if f0 == exon[1]:
				f0exons.append(exon)
			
			if f1 == exon[0]:
				f1exons.append(exon)
	
	best_f0_ac = None
	if f0exons != []:
		ac = Counter([exon[0] for exon in f0exons])
		accs = [a for a in ac]
		accc = [ac[a] for a in ac]

		maxi = accc.index(max(accc))

		best_f0_ac = accs[maxi]
	
	best_f1_do = None
	if f1exons != []:
		dc = Counter([exon[1] for exon in f1exons])
		dons = [a for a in dc]
		donc = [dc[a] for a in dc]

		maxi = donc.index(max(donc))

		best_f1_do = dons[maxi]
	

	return ((best_f0_ac, f0), (f1, best_f1_do))
	

def proc_snlist_sslist(strand, splc):
	## snlist is the ordered list of coordinates of successive splice sites involved in the event; 
	## sslist is the list of symbols (^, -) associated with each splice site in snlist, ordered according to snlist; 
	### Process the isoforms in the splice chai to collect coordinates and their corresponding symbols in snlist and sslist respec.
	snlist = []
	sslist = []
	for sc in splc.split(","):
		if sc == "":
			continue
		snlist.extend(re.sub(r'[\^\-]', ',', sc).split(',')[:-1])
		sslist.extend(re.sub(r'[\d]+', ',', sc).split(',')[1:])
	snlist = snlist
	snlist = [int(e) for e in snlist]

	sslist = sslist
	if strand == "+":
		sslist = [x for _, x in sorted(zip(snlist,sslist), key=lambda pair: pair[0])] ## sort the symbols according to the splice site coordinates
		snlist = sorted(snlist) ## sort the splice site coordinates now

	if strand == "-":
		sslist = [x for _, x in sorted(zip(snlist,sslist), reverse = True, key=lambda pair: pair[0])] ## sort the symbols according to the splice site coordinates
		snlist = sorted(snlist, reverse = True) ## sort the splice site coordinates now
	
	## combine snlist and sslist to get the master chain
	comb = []
	for (x,y) in zip(snlist, sslist): 
		comb.append(str(x)+y)
	
	## get the unique sequence of splice sites ## cannot use set() because it looses the ordered nature  
	combq = []
	for c in comb:
		if c not in combq:
			combq.append(c)

	scombq = "".join(combq) ## splice sequence as a string
	
	
	return scombq


def make_supp_tx_list_for_exp_isos(txs_supp, exp_iso):
	
	txs_supp = txs_supp.split(",")
	
	txs_supp_all = []
	for txs in txs_supp:
		txs_supp_all.extend(txs.split("/"))

	txs_supp_new = []
	for iso in exp_iso.split(","):
		idx = int(iso.split(".")[-1]) - 1
		txs_supp_new.append(txs_supp[idx])
	
	return (txs_supp_all, txs_supp_new)


def get_master_chain(strand, splc, txs_supp):
	
	## txs_supp --> list of supporting txs for each isoform expressed for an event -- for example [tx1/tx2/tx3, tx4/tx5, tx6] 

	pa5, pa3 = re.compile(r'(\d+\^){2,}'), re.compile(r'(\d+\-){2,}')
	scombq = proc_snlist_sslist(strand, splc)  										
	
	nscombq = scombq

	a5ss = [m.group() for m in re.finditer(pa5, scombq)]   ### find all instances of alt 5p ss and count the num of alt 5 ss involved in each such case
	#print("a5ss", a5ss)
	if a5ss != []:
		for g in a5ss:
			#print(g)
			ss, sg = [], []
			for s in g.split("^")[:-1]:
				for i in range(len(splc.split(","))):
					if s in splc.split(",")[i]:
						ss.append(s)
						sg.append(txs_supp[i].count("/") + 1)
			
			#print(ss, sg)
			maxi = sg.index(max(sg))
			best_s = ss[maxi]
			
			nscombq = nscombq.replace(g, best_s+"^")

	a3ss = [m.group() for m in re.finditer(pa3, scombq)]   ### find all instances of alt 5p ss and count the num of alt 5 ss involved in each such case
	#print("a3ss", a3ss)
	if a3ss != []:
		for g in a3ss:
			#print(g)
			ss, sg = [], []
			for s in g.split("-")[:-1]:
				for i in range(len(splc.split(","))):
					if s in splc.split(",")[i]:
						ss.append(s)
						sg.append(txs_supp[i].count("/") + 1)
			
			#print(ss, sg)
			maxi = sg.index(max(sg))
			best_s = ss[maxi]
			
			nscombq = nscombq.replace(g, best_s+"-")
	
	return nscombq

def get_asta_line(eid, astaf):
	
	lno = 0
	for e in eid.split("."):
		if "event" in e:
			lno = int(e.replace("event", "")) - 1

	line = ""
	with open(astaf, "r") as astah:
		line = astah.readlines()[lno]
	
	return line


def add_flanking_exons(f0exon, f1exon, master, strand):

	if strand == "+":
		if f0exon[0] != None:
			master = str(f0exon[0]) + "-" + str(f0exon[1]) + "^" + master
		if f0exon[0] == None:
			master = str(f0exon[1]) + "-" + master 
			
		if f1exon[1] != None:
			master = master + str(f1exon[0]) + "-" + str(f1exon[1]) + "^"
		if f1exon[1] == None:
			master = master + str(f1exon[0]) + "^"


	if strand == "-":
		if f0exon[0] != None:
			master = master + str(f0exon[1]) + "-" + str(f0exon[0]) + "^"
		if f0exon[0] == None:
			master = master + str(f0exon[1]) + "^" 
		
		if f1exon[1] != None:
			master = str(f1exon[1]) + "-" + str(f1exon[0]) + "^" + master 
		if f1exon[1] == None:
			master = str(f1exon[0]) + "-" + master

	return master

"""
- strand: 53331340-53331119^53327752-53327732^53327183-53327059^
+ strand: 24133943-24134054^24135746-24135875^
"""

def add_flanking_exons_by_tx(f0exon, f1exon, splc, strand):

	new_splc = []
	for sc in splc.split(","):
		if strand == "+":
			if f0exon[0] != None:
				sc = str(f0exon[0]) + "-" + str(f0exon[1]) + "^" + sc
			if f0exon[0] == None:
				sc = str(f0exon[1]) + "-" + sc 
				
			if f1exon[1] != None:
				sc = sc + str(f1exon[0]) + "-" + str(f1exon[1]) + "^"
			if f1exon[1] == None:
				sc = sc + str(f1exon[0]) + "^"


		if strand == "-":
			if f0exon[0] != None:
				sc = sc + str(f0exon[1]) + "-" + str(f0exon[0]) + "^"
			if f0exon[0] == None:
				sc = sc + str(f0exon[1]) + "^" 
			
			if f1exon[1] != None:
				sc = str(f1exon[1]) + "-" + str(f1exon[0]) + "^" + sc 
			if f1exon[1] == None:
				sc = str(f1exon[0]) + "-" + sc
		
		new_splc.append(sc)
	
	new_splc = ",".join(new_splc)

	return new_splc 

		
def get_effectv_lens(new_splc, readL, ancL):
	
	pinr = re.compile(r'\d+\^\d+\-')
	eff_lens = []

	for sc in new_splc.split(","):
		numJn = len(re.findall(pinr, sc))
		eff_lens.append((readL-ancL)*numJn)
	
	
	return eff_lens


def get_exons_and_regions(strand, master, lenT):

	pex = re.compile(r'\d+\-\d+\^')
	exons = []
	for m in re.finditer(pex, master):
		exon = m.group()
		exon = sorted([int(e) for e in exon.replace("^", "").split("-")])
		exons.append(exon)
	
	exon5, exon3 = exons[0], exons[-1]
	if strand == "-":
		exon3, exon5 = exons[0], exons[-1]

	return (exon5, exon3, exons[1:-1])


def get_introns_and_regions(strand, master, lenT):
	
	pinr = re.compile(r'\d+\^\d+\-')
	introns = []
	for m in re.finditer(pinr, master):
		inr = m.group()
		inr = sorted([int(e) for e in inr.replace("-", "").split("^")])
		introns.append(inr)

	return introns






















