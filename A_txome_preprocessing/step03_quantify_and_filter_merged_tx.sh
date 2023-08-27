module load salmon/1.4.0-gompi-2021a
module load python/3.10.4-gcccore-11.3.0

# make Salmon Index
MER="/scratch/user/uqktiwar/NCCSOP/merged"
IDX="/QRISdata/Q3294/NCCidx"
#mkdir $IDX
#mkdir $IDX"/allsamp"
#mkdir $IDX"/no07"

#salmon index -t $MER"/soxPN.f0.01.gtf.fa" -i $IDX"/allsamp" -p 36
#salmon index -t $MER"/soxPN.no07.f0.01.gtf.fa" -i $IDX"/no07" -p 36

# Run quantification with Salmon
FQ="/QRISdata/Q3294/Kanu/neural_crest/Fastq/merged"

SAL="/scratch/user/uqktiwar/NCCSALQ"
#mkdir $SAL


for n in 24707,ATCACG #24701,GTGGCC 24703,CGTACG 24704,GAGTGG 24705,ACTGAT 24706,ATTCCT 24707,ATCACG 24708,CGATGT 24709,TTAGGC 24710,TGACCA 24711,ACAGTG 24712,GCCAAT
 do
  IFS=","
  set $n
  mkdir $SAL"/"$1
  salmon quant -i $IDX"/allsamp" -l A \
  -1 $FQ"/"$1"_"$2"_R1.fastq.gz" \
  -2 $FQ"/"$1"_"$2"_R2.fastq.gz" \
  -p 36 --validateMappings --rangeFactorizationBins 4 --gcBias --seqBias -o $SAL"/"$1
 done

# Make a tx2gene file 
python createTx2gene.py $MER"/soxPN.f0.01.gtf" $MER"/soxPN.f0.01.tx2g.txt" stie
 
# get transcripts expressed according to user set TPM and samples thresholds
python filter_transcript_by_TPM.py -sqdir $SAL \
	-tx2g $MER"/soxPN.f0.01.tx2g.txt" \
	-p $MER"/samples_info.txt" \
	-tpm 0.5 -n 5 \
	-o $MER"/soxPN.tx2keep.txt"

# filter thddddddd'edlefe transcripts gtf and regenerate tx2gene file
python filter_transcripts_gtf.py $MER"/soxPN.tx2keep.txt" $MER"/soxPN.f0.01.gtf"
python createTx2gene.py $MER"/soxPN.f0.01.filt.gtf" $MER"/soxPN.f0.01.filt.tx2g.txt" stie

wc -l $MER"/soxPN.f0.01.filt.tx2g.txt"
wc -l $MER"/soxPN.tx2keep.txt"

