import sys
import argparse as sp
import re 
import subprocess 
import time

"""
Arguments 
arg1 - transcripts TPM table obtained after filtering format: TID\tSample1\tSample2 ... 
arg2 - transcripts GTF file obtained from stringtie 
"""

# get acceptable transcript ids 
valid_txids = []
lno = 0
with open(sys.argv[1], "r") as inph:
    for line in inph:
        lno += 1
        
        if lno == 1:
            # skip header
            continue

        txid = line.strip().split(",")[0]

        valid_txids.append(txid)


# read transcript gtf file, get txid of record and write the record if acceptable 
s=time.time()
with open(sys.argv[2].replace(".gtf", ".filt.gtf"), "w") as outh:
    with open(sys.argv[2], "r") as inph:
        for line in inph:
            if line[0] == "#":
                continue
            
            rec8 = line.strip().split("\t")[8]

            txid = rec8.replace('"', '').replace(";", "").split(" ")[3]

            if txid in valid_txids:
                outh.write(line)
e=time.time()
print(e-s)

"""
s = time.time()
with open(sys.argv[2].replace(".gtf", ".grep_filt.gtf"), "w") as outh:
    for txid in valid_txids:
        command = ["grep", "-wF", txid, sys.argv[2]]

        try:
            output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
            #print(output)
            outh.write(output)
        except subprocess.CalledProcessError as e:
            print("Error:", e.output)

e=time.time()
print(e-s)
"""
