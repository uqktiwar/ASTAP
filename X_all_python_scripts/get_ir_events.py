import sys
import re

## Input is asta OP file in GTF format
## chr1	Undefined	as_event	321032	322038	.	+	.	transcript_id "ENST00000599771,ENST00000601486"; gene_id "chr1:317730-322203W"; flanks "321032-,322038-"; structure "1^,2^"; splice_chain "321264^,321290^"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2"; 

## chr10	Undefined	as_event	73581633	73585650	.	-	.	transcript_id "MSTRG.3714.8,ENST00000394934,ENST00000394936/MSTRG.3714.12/MSTRG.3714.3/MSTRG.3714.4"; gene_id "chr10:73576055-73611126C"; flanks "73585650-,73581633^"; structure "0,1^2-,1^3-"; splice_chain ",73585594^73581770-,73585594^73581764-"; sources "Undefined,Undefined,Undefined"; NMorMRNA "null"; degree "9"; dimension "3_3"; 


outf1 = ".".join(sys.argv[1].split(".")[:-1]) + ".IR"
outf1b = ".".join(sys.argv[1].split(".")[:-1]) + ".IR.gtf"
outh1b = open(outf1b, "w")

outf2 = ".".join(sys.argv[1].split(".")[:-1]) + ".NotIR"
outf2b = ".".join(sys.argv[1].split(".")[:-1]) + ".NotIR.gtf"
outh2b = open(outf2b, "w")

n = 0
with open(outf1, "w") as outh1:
	with open(outf2, "w") as outh2:
		with open(sys.argv[1]) as inph:
			for line in inph:
				n += 1
				patt = line.split("\t")[8].split(" ")[7].replace('"', '').replace(';', '')
				flanks = line.split("\t")[8].split(" ")[5].replace('"', '').replace(';', '')
				flanks = re.sub('\d', '', flanks).replace("[", '-').replace("]", '^').split(",")
				
				if patt == '0,1^2-':
					#outh.write("event{}\n".format(n))
					outh1.write("event"+str(n)+".\n")
					outh1b.write(line)
					continue

				patt = patt.split(",")   ## [0, 1^2-, 1^3-] [1-2^4-8^9-11^14-16^17-, 3-8^9-10^15-16^17-,4-5^12-16^17-, 4-6^13-16^17-, 4-7^18-, 4-8^9-11^14-16^17-] 
				patt = ["-1"+flanks[0]+p+"1000"+flanks[1] for p in patt] ## patt = ["-1-01000^", "-1-1^2-1000^", "-1-1^3-1000^"] ["-1^1-2^4-8^9-11^14-16^17-1000^", "-1^3-8^9-10^15-16^17-1000^", "-1^4-5^12-16^17-1000^", "-1^4-6^13-16^17-1000^", "-1^4-7^18-1000^", "-1^4-8^9-11^14-16^17-1000^"] 
				
				## re.findall('\-1\-\d+\^|\d+\-\d+\^', '-1-1^2-3^100-')
				exons = [re.findall('^\-1\-\d+\^|\d+\-\d+\^', p) for p in patt] 
				introns = [re.findall('^\-1\^\d+\-|\d+\^\d+\-', p) for p in patt]
				
				if [] in introns:
					outh1.write("event"+str(n)+".\n")
					outh1b.write(line)
					continue
				
				foundir = 0
				for i in range(len(exons)):
					for j in range(len(exons)):
						if i == j:
							continue
						
						for inr in introns[i]:
							inrc = [int(s) for s in inr[:-1].split("^")]
							for exn in exons[j]:
								exnc = [0,0]
								if len(exn[:-1].split("-")) == 3:
									exnc = [-1, int(exn[:-1].split("-")[2])]
								else:
									exnc = [int(s) for s in exn[:-1].split("-")]
								
								if inrc[0] > exnc[0] and inrc[1] < exnc[1]:
									#outh.write("event{}\n".format(n))
									foundir = 1
									outh1.write("event"+str(n)+".\n")
									outh1b.write(line)
									break
							
							if foundir:
								break
					
						if foundir:
							break
					
					if foundir:
						break
				
				if foundir == 0:
					outh2.write("event"+str(n)+".\n")
					outh2b.write(line)

outh1b.close()
outh2b.close()
"""
flanks "6856645^,6835982^"; structure "1-2^4-8^9-11^14-16^17-,3-8^9-10^15-16^17-,4-5^12-16^17-,4-6^13-16^17-,4-7^18-,4-8^9-11^14-16^17-"
f0---1####2---------4#####################################################################8--------9####################11-------------------------14#################16-------------17##############f1
f0---------------3########################################################################8--------9###############10------------------------------------15###########16-------------17##############f1
f0------------------4##################5----------------------------------------------------------------------------------------12####################################16-------------17##############f1
f0------------------4###############################6--------------------------------------------------------------------------------------13#########################16-------------17##############f1
f0------------------4################################################7----------------------------------------------------------------------------------------------------------------------18#######f1
f0------------------4#####################################################################8--------9####################11-------------------------14#################16-------------17##############f1



### 0*1 ----- 2*3 ----------------------------- 8,
	0*1 ----- 2*3 ---- 4*5 ---- 6*7 ----------- 8,
	0*1 ----- 2*3 ---- 4**********7 ----------- 8,
	0*1 ----- 2*3 --------------6*7 ----------- 8,
	0*1 ----- 2**********5 ---- 6*7 ----------- 8,
	0***********3 ---- 4*5 ---- 6*7 ----------- *						
"""

				
				
			



