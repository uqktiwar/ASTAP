module load java/11.0.18
module load htslib/1.15.1-gcc-11.3.0
module load samtools/1.13-gcc-10.3.0
module load python/3.10.4-gcccore-11.3.0-bare
module load pysam/0.16.0.1-gcc-10.3.0

BASE="/scratch/user/uqktiwar/NCCSPLICE"
MER="/scratch/user/uqktiwar/NCCSOP/merged"
ASTA="/home/uqktiwar/ASTAP"

## 1.) Make the NotIR req events
mkdir $BASE"/NOTIR_req"
for n in {1..3}
   do
        python $ASTA"/asta_to_bed.requant.py" $BASE"/soxPN.NotIR.allD.events2requant.modFlanks.out" $BASE"/NOTIR_req/soxPN.NotIR."$n
   done

## 2.) Make the IR req events
mkdir $BASE"/IR_req"
python $ASTA"/asta_to_bed.requant.py" $BASE"/soxPN.IR.allD.events2requant.modFlanks.out" $BASE"/IR_req/soxPN.IR"

## 3.) Copy the exon database
cp $BASE"/IR/asg.exon_database.bed.gz" $BASE"/IR_req/asg.exon_database.bed.gz"
tabix -p bed $BASE"/IR_req/asg.exon_database.bed.gz"


## 4.) Pre-prop the IR events
python $ASTA"/get_ir_event_features.py" $BASE"/IR_req" "soxPN.IR.bed"
python $ASTA"/find_intron_exon_ovlaps.py" $BASE"/IR_req" "all_retained_introns.bed" "ret_introns_regions.out" 0

python $ASTA"/get_all_junctions.py" $BASE"/IR_req" "soxPN.IR.bed" "soxPN"

for n in {1..3}
  do
      cp $BASE"/IR_req/ret_introns_regions.out" $BASE"/IR_req/ret_introns_regions."$n".out"
      cp $BASE"/IR_req/soxPN.IR.junctions.bed" $BASE"/IR_req/soxPN.IR.junctions."$n".bed"

  done




