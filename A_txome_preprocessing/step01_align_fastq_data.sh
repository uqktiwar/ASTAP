#!/bin/bash 

# Add STAR executable to you $PATH
module load star/2.7.9a-gcc-10.3.0
module load samtools

# Specify directory containing genome fasta file(s). Check STAR manual to see which fasta file(s) you need for genomeGenerate step
GENFA="/QRISdata/Q3294/Kanu/30days_flite/genomes/human/hg38"

# Specify path to store STAR index, make sure the folder exists but is empty
GEN="/QRISdata/Q3294/Kanu/STAR_indexes/hg38_99"
#rm -rf $GEN 
#mkdir $GEN

# Path to store aligned BAM output files
AOP="/scratch/user/uqktiwar/NCCAOP"
mkdir $AOP

# Path to input FASTQ files
FQ="/QRISdata/Q3294/Kanu/neural_crest/Fastq/merged"

# Generate STAR Index. The sjdbOverhang param = fastq read length - 1  # Skip this step if you already have an index file

#STAR --runThreadN 32 \
#--runMode genomeGenerate \
#--genomeDir $GEN \
#--genomeFastaFiles $GENFA"/hg38.primary_assembly.fa" \
#--sjdbGTFfile $GENFA"/Homo_sapiens.GRCh38.102.chr.gtf" \
#--sjdbOverhang 99


# Align Fastq files. The names of the fatqs are stored in fastq_names_all.txt file. filenames are formmatted as abc_deg_R1.fastq.gz  
# Only change runThreadN  

IFS="_"
while read -ra line
	do 
		echo "${line[0]}_${line[1]}_R1.fastq.gz" 
		echo "${line[0]}_${line[1]}_R2.fastq.gz" 

		STAR --runMode alignReads --runThreadN 32 --genomeDir "$GEN" \
		--genomeLoad NoSharedMemory --readFilesIn "$FQ/${line[0]}_${line[1]}_R1.fastq.gz" "$FQ/${line[0]}_${line[1]}_R2.fastq.gz" \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix "$AOP/${line[0]}" \
		--outSAMattributes NH HI NM MD AS nM jM jI XS \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 37125051484

		samtools index -b -@ 32 "$AOP/${line[0]}Aligned.sortedByCoord.out.bam"

	done < /home/uqktiwar/ASTAP/fastq_names_all.txt



