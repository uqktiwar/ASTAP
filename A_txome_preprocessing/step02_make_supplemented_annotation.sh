#!/bin/bash

# 1.) Run stringtie on each sample 

# Specify directory containing aligned BAM files
AOP=/scratch/user/uqktiwar/NCCAOP

# SPecify output directory to store reconstructed transcripts by stringtie
SOP=/scratch/user/uqktiwar/NCCSOP
mkdir $SOP

# Specify path to stringtie directory containing stringtie executable
STIE="/home/uqktiwar/ASTAP/stringtie"

# Add stringtie path to your $PATH
export PATH=$PATH:$STIE

# specify reference file -- same file used for genome generate, mapping with STAR

REF="/QRISdata/Q3294/Kanu/30days_flite/genomes/human/hg38/Homo_sapiens.GRCh38.102.chr.gtf"

# reconstruct tx for each sample
for n in {24703..24712}
  do
	 
	 stringtie --rf $AOP/$n"Aligned.sortedByCoord.out.bam" -o $SOP/$n/$n".f0.01.gtf" -l $n -p 36 -G $REF -f 0.01 

  done

stringtie --rf $AOP"/24701Aligned.sortedByCoord.out.bam" -o $SOP"/24701/24701.f0.01.gtf" -l "24701" -p 36 -G $REF -f 0.01

# Run stringtie-merge
mkdir $SOP/merged

#stringtie --merge -f 0.01 -G $REF -p 36 -o $SOP"/merged/soxPN.no07.f0.01.gtf" assembly_GTF_list_no07.txt
stringtie --merge -f 0.01 -G $REF -p 36 -o $SOP"/merged/soxPN.f0.01.gtf" assembly_GTF_list.txt

# Run salmon quant using the merged annotation and BAM files
# generate merged tx fasta
/home/uqktiwar/ASTAP/gffread/./gffread -g "/QRISdata/Q3294/Kanu/30days_flite/genomes/human/hg38/hg38.primary_assembly.fa" -w $SOP"/merged/soxPN.f0.01.gtf.fa" --gtf $SOP"/merged/soxPN.f0.01.gtf"
#/home/uqktiwar/ASTAP/gffread/./gffread -g "/QRISdata/Q3294/Kanu/30days_flite/genomes/human/hg38/hg38.primary_assembly.fa" -w $SOP"/merged/soxPN.no07.f0.01.gtf.fa" --gtf $SOP"/merged/soxPN.no07.f0.01.gtf"



