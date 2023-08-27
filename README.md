# ASTA-P 
### A pipeline for the detection, quantification and statistical analysis of complex alternative splicing events
<div align="justify"> 
Alternative splicing (AS) is an important regulatory mechanism for modulating gene expression or protein function across different biological conditions. Vast amounts of RNA-seq data surveying different physiological and disease conditions are continuously being generated and analysed to advance our knowledge of AS mechanisms, and it is now estimated that more than 95% of human genes express alternative isoforms[1]. <br><br>Traditionally, splicing analyses using RNA-seq data involve pairwise comparisons between isoforms of a gene, uncovering a set of commonly analysed “classical” binary splicing patterns, namely: cassette exon skipping, alternative 5’ or 3’ splice site (SS) usage, mutually exclusive exons, and intron retention [2-4]. However, overlapping instances of these patterns are often observed, indicating that a more complex cumulative AS structure would be a better descriptor of the splicing variation observed for some genes. <br><br>We present ASTA-P, a pipeline for profiling, quantification, and differential splicing analysis of tissue-specific, arbitrarily complex alternative splicing patterns. We discover novel events by supplementing existing annotation with reconstructed transcripts and use spliced RNA-seq reads to quantify splicing changes accurately based on their unique assignments.
For each gene, arbitrarily complex splicing events are mined as bubbles from the alternative splicing graph using the ASTALAVISTA [5] algorithm, followed by quantification with a custom script, and modelling the events’ counts using the Dirichlet-multinomial regression framework implemented in the R-package, DRIMseq [6].
</div>

## Pipeline Schematic
<img align="center" src="https://github.com/uqktiwar/ASTAP/assets/10694707/83d342b2-110e-4fb6-a50d-7d178667212e">

## References
<div align="justify">
1.	Q. Pan, O. Shai, L. J. Lee, B. J. Frey, and B. J. Blencowe, “Deep surveying of alternative splicing complexity in the human transcriptome by high-throughput sequencing,” Nat. Genet., vol. 40, no. 12, pp. 1413–1415, 2008.<br>
2.	S. Shen et al., “rMATS: Robust and flexible detection of differential alternative splicing from replicate RNA-Seq data,” Proc. Natl. Acad. Sci., vol. 111, no. 51, pp. E5593--E5601, 2014.<br>
3.	J. L. Trincado et al., “SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions,” Genome Biol., vol. 19, no. 1, p. 40, 2018.<br>
4.	Katz, Y., Wang, E., Airoldi, E. et al. Analysis and design of RNA sequencing experiments for identifying isoform regulation. Nat Methods 7, 1009–1015 (2010).<br>
5.	S. Foissac and M. Sammeth, “ASTALAVISTA : dynamic and flexible analysis of alternative splicing events in custom gene datasets,” Nucleic Acids Res., vol. 35, pp. 297–299, 2007.<br>
6.	M. Nowicka and M. D. Robinson, “DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics,” F1000Research, vol. 5, p. 1356, 2016.<br>
</div>

