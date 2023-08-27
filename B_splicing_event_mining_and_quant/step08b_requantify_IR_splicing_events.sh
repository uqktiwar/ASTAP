module load python/3.10.4-gcccore-11.3.0
module load pysam/0.16.0.1-gcc-10.3.0 

# Quantify NotIR events' features
BAM="/scratch/user/uqktiwar/NCCAOP"
BASE="/home/uqktiwar/ASTAP"
SPLC="/scratch/user/uqktiwar/NCCSPLICE/IR_req"

OP=$SPLC"/IR1"
mkdir $OP
python $BASE/find_ovlappin_introns_and_quant_ir.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM \
	-bams 24701 24703 24704 24705 24706 24707 24708 24709 24710 24711 24712 \
	-s "Aligned.sortedByCoord.out" -l rf

python $BASE/find_skipping_jn_counts_for_IR.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM \
	-bams 24701 24703 24704 24705 24706 24707 24708 24709 24710 24711 24712 \
	-s "Aligned.sortedByCoord.out" -l rf

python $BASE/quant_junctions_for_IR.1.py -w $SPLC -i soxPN.IR.junctions.1.bed -o $OP -b $BAM \
	-bams 24701 24703 24704 24705 24706 24707 24708 24709 24710 24711 24712 \
	-s "Aligned.sortedByCoord.out" -l rf


<<comment
OP=$SPLC"/IR2"
mkdir $OP
python $BASE/find_ovlappin_introns_and_quant_ir.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM -bams 24705 24706 24707 24708 -s "Aligned.sortedByCoord.out" -l rf
python $BASE/find_skipping_jn_counts_for_IR.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM -bams 24705 24706 24707 24708 -s "Aligned.sortedByCoord.out" -l rf

python $BASE/quant_junctions_for_IR.1.py -w $SPLC -i soxPN.IR.junctions.1.bed -o $OP -b $BAM -bams 24705 24706 24707 24708 -s "Aligned.sortedByCoord.out" -l rf


OP=$SPLC"/IR3"
mkdir $OP
python $BASE/find_ovlappin_introns_and_quant_ir.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM -bams 24709 24710 24711 24712 -s "Aligned.sortedByCoord.out" -l rf
python $BASE/find_skipping_jn_counts_for_IR.1.py -w $SPLC -i ret_introns_regions.1.out -o $OP -b $BAM -bams 24709 24710 24711 24712 -s "Aligned.sortedByCoord.out" -l rf

python $BASE/quant_junctions_for_IR.1.py -w $SPLC -i soxPN.IR.junctions.1.bed -o $OP -b $BAM -bams 24709 24710 24711 24712 -s "Aligned.sortedByCoord.out" -l rf

comment



