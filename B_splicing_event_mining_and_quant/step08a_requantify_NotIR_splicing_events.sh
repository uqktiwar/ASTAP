module load python/3.10.4-gcccore-11.3.0
module load pysam/0.16.0.1-gcc-10.3.0 

# Quantify NotIR events' features
BAM="/scratch/user/uqktiwar/NCCAOP"
BASE="/home/uqktiwar/ASTAP"
SPLC="/scratch/user/uqktiwar/NCCSPLICE/NOTIR_req"

OP=$SPLC"/NIR1"
mkdir $OP
python $BASE"/countsADJ.NotIR.1.py" -w $SPLC -i "soxPN.NotIR.1.bed" -o $OP \
	-b $BAM -bams 24701 24703 24704 24705 24706 24707 24708 24709 24710 24711 24712 \
	-s "Aligned.sortedByCoord.out" -l rf

<<comment
OP=$SPLC"/NIR2"
mkdir $OP
python $BASE"/countsADJ.NotIR.1.py" -w $SPLC -i "soxPN.NotIR.1.bed" -o $OP \
        -b $BAM -bams 24705 24706 24707 24708 \
        -s "Aligned.sortedByCoord.out" -l rf

OP=$SPLC"/NIR3"
mkdir $OP
python $BASE"/countsADJ.NotIR.1.py" -w $SPLC -i "soxPN.NotIR.1.bed" -o $OP \
        -b $BAM -bams 24709 24710 24711 24712 \
        -s "Aligned.sortedByCoord.out" -l rf
comment


