library(dplyr)

setwd(wd <- "D:/NeuralCrest")
dge.dir <- file.path(wd, "salQ")
ds.dir <- file.path(wd, "splice/NotIR/drimOP")

dge.bpt <- "soxPN"
dge.col <- "log2FoldChange"

splc.bpt <- "SoxNP"
cells <- c("SOXN", "SOXP")

dge.bpt.df <- read.table(file.path(dge.dir, paste(dge.bpt, "DGE.out", sep = ".")), header = T, sep = "\t", stringsAsFactors = F)
dge.bpt.df <- dge.bpt.df[dge.bpt.df$padj < 0.05,]

##EID	CHR	STRAND	GID	EGID	HGNC	NFLANK.5P	NFLANK.3P	ISOS.EXP	NEW.SPLC	OD	RD	
##feature_id	gene	transcript	HSC	CMP	CLP	HSC.FC	CMP.FC	CLP.FC	HSC.Nr	CMP.Nr	CLP.Nr

splc.df <- read.table(file.path(ds.dir, paste(splc.bpt, "final.Med30.out", sep = ".")), header = T, sep = "\t", stringsAsFactors = F)
splc.df[, paste("delPsi", splc.bpt, sep = ".")] <- splc.df[, cells[2]] - splc.df[, cells[1]]
splc.df[, "sig.DGE"] <- 0
splc.df[, "lfc.DGE"] <- 0
  
dge.sub.idx <- which(dge.bpt.df[ ,"padj"]<0.05 & abs(dge.bpt.df[ , dge.col]) >= 1.5)
dge.sub.hgnc <- dge.bpt.df$HGNC[dge.sub.idx]
splc.df[splc.df$HGNC %in% dge.sub.hgnc, "sig.DGE"] <- 1
  
###print(which(splc.df$HGNC %in% dge.bpt.df$HGNC))
for(i in dge.sub.idx){splc.df[which(splc.df$HGNC == dge.bpt.df$HGNC[i]), "lfc.DGE"] <- dge.bpt.df[i, dge.col]}
  
write.table(splc.df, file.path(ds.dir, paste(splc.bpt, "final.Med30.w.dge.out", sep = ".")), quote = F, row.names = F, sep = "\t")





