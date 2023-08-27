module load python/3.10.4-gcccore-11.3.0

BASE="/scratch/user/uqktiwar/NCCSPLICE"
MER="/scratch/user/uqktiwar/NCCSOP/merged"
ASTA="/home/uqktiwar/ASTAP"

## filter isoforms within the NotIR events ## REMEMBER TO CHANGE THE FILENAMES OF THE NOTIR INPUT COUNT FILES TO MATCH THE CFPJN OPTION
<<comment1
for n in {1..3}
  do
        mv $BASE"/NOTIR/NIR"$n"/soxPN.NotIR.1.bed.feat.countsAdj.out" $BASE"/NOTIR/NIR"$n"/soxPN.NotIR.feat."$n".countsAdj.out"
        mv $BASE"/NOTIR/NIR"$n"/soxPN.NotIR.1.bed.countsAdj.out" $BASE"/NOTIR/NIR"$n"/soxPN.NotIR.iso."$n".countsAdj.out"
  done

python $ASTA"/filter_isoforms.py" -w $BASE -mjd $BASE/NOTIR \
                          -pdf pdata.nir.txt -cfpJn soxPN.NotIR.feat \
                          -iff NOTIR/soxPN.NotIR.1.bed -etype NotIR \
			  -ft 5 -mins 5 \
			  -op $BASE/soxPN.NotIR


## filter isoforms within IR events
python $ASTA"/filter_isoforms.py" -w $BASE -mjd $BASE/IR \
                          -pdf pdata.ir.txt -cfpJn soxPN.IR.junctions \
                          -cfpRi ret_introns_regions \
                          -iff IR/all_ir_events_features.out \
                          -etype IR -op $BASE/soxPN.IR \
			  -ft 5 -mins 5


## make eevent to gene mapping files for all events (all dimensions) use tx2g file already created in Step02
python $ASTA"/event2gene.py" $BASE"/soxPN.asta.NotIR.gtf" $MER"/soxPN.f0.01.tx2g.txt" $BASE"/soxPN.NotIR.allD.e2g.out"
python $ASTA"/event2gene.py" $BASE"/soxPN.asta.IR.gtf" $MER"/soxPN.f0.01.tx2g.txt" $BASE"/soxPN.IR.allD.e2g.out"

## remove low isoforms and make new splice chains
python $ASTA"/rm_low_isoforms_and_redo_splc.py" $BASE"/soxPN.NotIR.allD.e2g.out" $BASE"/soxPN.NotIR.FeatExprsn.out" $BASE"/soxPN.NotIR.newSplc.out"
python $ASTA"/rm_low_isoforms_and_redo_splc.py" $BASE"/soxPN.IR.allD.e2g.out" $BASE"/soxPN.IR.FeatExprsn.out" $BASE"/soxPN.IR.newSplc.out"
comment1

## remove overlapping events and mark events which need to be requantified on account of changing flanks
python $ASTA"/rm_ovlapping_events_and_mark_for_requant.py" $BASE"/soxPN.NotIR.newSplc.out" $BASE"/soxPN.NotIR.allD"
python $ASTA"/rm_ovlapping_events_and_mark_for_requant.py" $BASE"/soxPN.IR.newSplc.out" $BASE"/soxPN.IR.allD"


