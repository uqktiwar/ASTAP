import pandas as pd
import os
import sys
import argparse as ap

parser = ap.ArgumentParser(
                    prog='filter_transcripts.py',
                    description='This script filters the reconstructed transcripts according to TPM thresholds set by user.',)

parser.add_argument('-sqdir', '--salmon_quant_dir', help = "Base directory where salmon quant sub-directories for each RNA_seq sample are stored")
parser.add_argument('-tx2g', '--tx2g_file', help = "Full path to tx2gene mapping file")
parser.add_argument('-p', '--samples_file', help = "Full path to samples info file. This is a tsv file containing atleast a Sample and Group column.")
parser.add_argument('-tpm', '--tpm_threshold', type = float, help = "TPM threhold to use for filtering transcripts")
parser.add_argument('-n', '--min_samples', type = int, help = "min samples in which a tx should be expressed at or above TPM threhold")
parser.add_argument('-o', '--op_file', help = "Full path to output file name to store filtered transcripts")

args = parser.parse_args()

def read_salmon_quant_files(quant_dir, samples, tx2g):
    """
    quant_dir: base directory where salmon quant files are stored
    samples: sample subdirectories' names where salmon quant.sf files are stored
    tx2g: tx2gene database of transcripts
    """
    
    for s in samples:
        osamp = pd.read_csv(os.path.join(quant_dir, str(s), "quant.sf"), sep = "\t")
        osamp.columns = ["TID", "Length", "EffectiveLength", s, "NumReads"]
        osamp = osamp.loc[:, ["TID", s]]
        
        tx2g = pd.merge(tx2g, osamp, how="left", left_on="TID", right_on="TID")
    
    return tx2g.loc[:, ["TID"]+samples]


tx2g = pd.read_csv(args.tx2g_file, sep="\t")

samples = pd.read_csv(args.samples_file, sep = "\t")

final = read_salmon_quant_files(quant_dir=args.salmon_quant_dir, samples = samples.Sample.to_list(), tx2g=tx2g)

final.set_index("TID", inplace=True)

final_sub = final.loc[(final >= args.tpm_threshold).sum(axis = 1) >= args.min_samples, :]

print("Keeping", final_sub.shape[0], "out of", final.shape[0], "transcripts with atleast", args.tpm_threshold, "TPM in", args.min_samples, "samples....")

final_sub.loc[:, "TID"] = final_sub.index.to_list()
final_sub.loc[:, ["TID"] + samples.Sample.to_list()].to_csv(args.op_file, index=False)



