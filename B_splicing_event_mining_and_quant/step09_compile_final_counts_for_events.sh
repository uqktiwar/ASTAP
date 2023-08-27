module load java/11.0.18
module load htslib/1.15.1-gcc-11.3.0
module load samtools/1.13-gcc-10.3.0
module load python/3.10.4-gcccore-11.3.0-bare
module load pysam/0.16.0.1-gcc-10.3.0

BASE="/scratch/user/uqktiwar/NCCSPLICE"
MER="/scratch/user/uqktiwar/NCCSOP/merged"
ASTA="/home/uqktiwar/ASTAP"


## For NotIR events
#usage: compile_iso_level_counts.nir.py [-h] [-w W] [-mjd1 MASTERDIR1]
#                                       [-mjd2 MASTERDIR2] [-e2gf E2GFILE]
#                                       [-csubd COUNTSUBDIRS [COUNTSUBDIRS ...]]
#                                       [-cfp1 JNCOUNTFILEPREFIX1]
#                                       [-cfp2 JNCOUNTFILEPREFIX2]
#                                       [-op OUTPREFIX]
<<comment1
for n in 1
  do

        mv $BASE"/NOTIR_req/NIR"$n"/soxPN.NotIR."$n".bed.feat.countsAdj.out" $BASE"/NOTIR_req/NIR"$n"/soxPN.NotIR.feat."$n".countsAdj.out"
        mv $BASE"/NOTIR_req/NIR"$n"/soxPN.NotIR."$n".bed.countsAdj.out" $BASE"/NOTIR_req/NIR"$n"/soxPN.NotIR.iso."$n".countsAdj.out"

  done


python $ASTA"/compile_iso_level_counts.byIso.py" -w $BASE -e2gf soxPN.NotIR.allD.final.e2g.out \
-mjd1 NOTIR -mjd2 NOTIR_req \
-csubd NIR1 NIR2 NIR3 \
-cfp1 soxPN.NotIR.iso -cfp2 soxPN.NotIR.iso -op soxPN.NotIR.final

python $ASTA"/compile_iso_level_counts.byIso.py" -w $BASE -e2gf soxPN.NotIR.allD.final.e2g.out -mjd1 NOTIR -csubd NIR1 NIR2 NIR3 -cfp1 soxPN.NotIR.iso -op soxPN.NotIR.orig
comment1

#usage: compile_iso_level_counts.byFeat.py [-h] [-w W] [-mjd1 MASTERDIR1]
#                                          [-mjd2 MASTERDIR2]
#                                          [-pdf1 PDATAFILE1]
#                                          [-pdf2 PDATAFILE2]
#                                          [-cfpJn1 JNCOUNTFILEPREFIX1]
#                                          [-cfpJn2 JNCOUNTFILEPREFIX2]
#                                          [-cfpRi1 [RICOUNTFILEPREFIX1]]
#                                          [-cfpRi2 [RICOUNTFILEPREFIX2]]
#                                          [-iff INPUTFEATUREFILE]
#                                          [-e2gf E2GFILE] [-etype ETYPE]
#                                          [-op OUTPREFIX]


## only for IR events here 
python $ASTA"/compile_iso_level_counts.byFeat.py" -w $BASE -mjd1 IR -mjd2 IR_req\
                                                  -pdf1 $BASE/pdata.ir.txt -pdf2 $BASE/pdata.irr.txt\
                                                  -cfpJn1 soxPN.IR.junctions -cfpJn2 soxPN.IR.junctions \
                                                  -cfpRi ret_introns_regions -iff1 IR/all_ir_events_features.out\
                                                  -iff2 IR_req/all_ir_events_features.out -e2gf soxPN.IR.allD.final.e2g.out\
                                                  -etype IR -op soxPN.IR.final

python $ASTA"/compile_iso_level_counts.byFeat.py" -w $BASE -mjd1 IR -pdf1 $BASE/pdata.ir.txt -cfpJn1 soxPN.IR.junctions\
                                                  -cfpRi ret_introns_regions -iff1 IR/all_ir_events_features.out\
                                                  -e2gf soxPN.IR.allD.final.e2g.out -etype IR -op soxPN.IR.orig




