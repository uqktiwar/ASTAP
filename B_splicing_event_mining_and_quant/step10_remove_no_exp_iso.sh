module load java/11.0.18
module load htslib/1.15.1-gcc-11.3.0
module load samtools/1.13-gcc-10.3.0
module load python/3.10.4-gcccore-11.3.0-bare
module load pysam/0.16.0.1-gcc-10.3.0

BASE="/scratch/user/uqktiwar/NCCSPLICE"
MER="/scratch/user/uqktiwar/NCCSOP/merged"
ASTA="/home/uqktiwar/ASTAP"


python $ASTA"/correct_zero_exp_iso_counts.py" $BASE soxPN.NotIR.FeatExprsn.out soxPN.NotIR.final.IsoExprsn.out
mv $BASE/soxPN.NotIR.final.IsoExprsn.out.zero $BASE/soxPN.NotIR.final.IsoExprsn.out

python $ASTA"/correct_zero_exp_iso_counts.py" $BASE soxPN.NotIR.FeatExprsn.out soxPN.NotIR.orig.IsoExprsn.out
mv $BASE/soxPN.NotIR.orig.IsoExprsn.out.zero $BASE/soxPN.NotIR.orig.IsoExprsn.out

python $ASTA"/correct_zero_exp_iso_counts.py" $BASE soxPN.IR.FeatExprsn.out soxPN.IR.final.IsoExprsn.out
mv $BASE/soxPN.IR.final.IsoExprsn.out.zero $BASE/soxPN.IR.final.IsoExprsn.out

python $ASTA"/correct_zero_exp_iso_counts.py" $BASE soxPN.IR.FeatExprsn.out soxPN.IR.orig.IsoExprsn.out
mv $BASE/soxPN.IR.orig.IsoExprsn.out.zero $BASE/soxPN.IR.orig.IsoExprsn.out


