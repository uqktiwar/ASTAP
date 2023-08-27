module load htslib/1.12-gcc-10.3.0
module load samtools/1.13-gcc-10.3.0
module load java/11.0.18
module load python/3.10.4-gcccore-11.3.0
module load pysam/0.16.0.1-gcc-10.3.0

MER="/scratch/user/uqktiwar/NCCSOP/merged"
AOP="/scratch/user/uqktiwar/NCCSPLICE"
mkdir $AOP

# 1. Find Splicing events from filtered GTF
ASTA="/home/uqktiwar/ASTAP/astalavista-4.0.1/bin"

$ASTA/astalavista -t asta -i $MER"/soxPN.f0.01.filt.gtf" -e [ASI] -d 0 -o $AOP"/soxPN.asta.gtf.gz"
gunzip $AOP"/soxPN.asta.gtf.gz"

# 2. Separate events into IR and NOTIR events
mkdir $AOP"/NOTIR"
mkdir $AOP"/IR"

python get_ir_events.py $AOP"/soxPN.asta.gtf"

# 3. Convert events into bed format
python asta_to_bed.py $AOP"/soxPN.asta.NotIR.gtf" $AOP"/NOTIR/soxPN.NotIR"
mv $AOP"/NOTIR/soxPN.NotIR.bed" $AOP"/NOTIR/soxPN.NotIR.1.bed"
for n in {2..7}
  do
     cp $AOP"/NOTIR/soxPN.NotIR.1.bed" $AOP"/NOTIR/soxPN.NotIR."$n".bed"
  done

python asta_to_bed.py $AOP"/soxPN.asta.IR.gtf" $AOP"/IR/soxPN.IR"


# 4. Make and exon database for use during IR event quantification
./gtf2bed.sh $MER"/soxPN.f0.01.gtf" > $MER"/soxPN.f0.01.bed"    ## use the full gtf file to make the exon database 
python get_all_exons_from_bed.py $MER"/soxPN.f0.01.bed" $AOP"/asg.exon_database.bed"

# 5. Sort and Index the exon database for fast access from within python
sort -k1,1 -k2,2n $AOP"/asg.exon_database.bed" > $AOP"/asg.exon_database.s.bed"
mv $AOP"/asg.exon_database.s.bed" $AOP"/asg.exon_database.bed"
cp $AOP"/asg.exon_database.bed" $AOP"/IR/asg.exon_database.bed"
bgzip $AOP"/IR/asg.exon_database.bed"
tabix -p bed $AOP"/IR/asg.exon_database.bed.gz"

# 6. Some additional Pre-processing of the IR events 
python get_ir_event_features.py $AOP"/IR" "soxPN.IR.bed"

python find_intron_exon_ovlaps.py $AOP"/IR" "all_retained_introns.bed" "ret_introns_regions.out" 0
mv $AOP"/IR/ret_introns_regions.out" $AOP"/IR/ret_introns_regions.1.out"

python get_all_junctions.py $AOP"/IR" "soxPN.IR.bed" "soxPN"
mv $AOP"/IR/soxPN.IR.junctions.bed" $AOP"/IR/soxPN.IR.junctions.1.bed"

for n in {2..7}
  do
     cp $AOP"/IR/ret_introns_regions.1.out" $AOP"/IR/ret_introns_regions."$n".out"
     cp $AOP"/IR/soxPN.IR.junctions.1.bed" $AOP"/IR/soxPN.IR.junctions."$n".bed"
  done



