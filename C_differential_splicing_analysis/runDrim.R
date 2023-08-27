setwd(wd <- "/scratch/user/uqktiwar/NCCSPLICE/NOTIR")

library(DRIMSeq)
library("dplyr", verbose = F)
library(tidyr)
library(stageR)
source("getNrCounts.R")

opDir <- file.path(wd, "drimOP")

mingrp <- 5

OP <- "SoxNP"
cells <- c("SOXN", "SOXP")
repr.Samps <- c("N24712", "P24701")
samples <- read.table(file.path(wd, "pdata.nir.txt"), header = TRUE, stringsAsFactors = F)
e2g <- read.table(file.path(wd, "soxPN.NotIR.allD.final.e2g.out"), header = T, sep = "\t", stringsAsFactors = F)

##### Read countfiles, rename event ids by dimension, rbind counts from different dimensions ###### 
txiT <- read.table(file.path(wd, "soxPN.NotIR.final.IsoExprsn.out"), sep = "\t", header = T, stringsAsFactors = F)
txiT <- txiT[, c("gene_id", "feature_id", c(samples$Sample))]

nrcountsDFs <- getNrCounts(txiT, samples, "gene_id", mingrp, T)
nrcMin <- nrcountsDFs[[1]]
colnames(nrcMin) <- c("EID", paste(unique(samples$CellType), "Nr", sep = "."))
saveRDS(nrcMin, file.path(opDir, paste(OP, "nrcMin", sep = ".")))
nrcMin <- readRDS(file.path(opDir, paste(OP, "nrcMin", sep = ".")))
keep <- nrcMin$EID[rowSums(nrcMin[,2:ncol(nrcMin)] >= 20)==2]

nrcounts <- nrcountsDFs[[2]]
colnames(nrcounts) <- c("EID", paste(unique(samples$CellType), "Nr", sep = "."))
saveRDS(nrcounts, file.path(opDir, paste(OP, "nrcounts", sep = ".")))
nrcounts <- readRDS(file.path(opDir, paste(OP, "nrcounts", sep = ".")))
  
txcountsDFs <- getNrCounts(txiT, samples, "feature_id", mingrp, T)
txcounts <- txcountsDFs[[2]]
colnames(txcounts) <- c("feature_id", paste(unique(samples$CellType), "FC", sep = "."))
saveRDS(txcounts, file.path(opDir, paste(OP, "txcounts", sep = ".")))
txcounts <- readRDS(file.path(opDir, paste(OP, "txcounts", sep = ".")))

############################################################################################################################
txiT <- txiT[txiT$gene_id %in% keep,]
samples.mod <- samples
colnames(samples.mod)[1:2] <- c("sample_id", "CT")
samples.mod$CT <- factor(samples$CellType, levels = cells)

drimObj <- dmDSdata(counts = txiT, samples = samples.mod)
drimObj <- dmFilter(drimObj, min_samps_feature_expr = 5,
                     min_feature_expr = 5)

model.full <- model.matrix(~CT , data = samples.mod)

drimObj <- dmPrecision(drimObj, design = model.full, verbose = 2, add_uniform = T)
drimObj <- dmFit(drimObj, design = model.full, verbose = 2, add_uniform = T)
saveRDS(drimObj, file.path(opDir, paste(OP, "dobj.Min20.MF5.fit", sep = ".")))

props <- proportions(drimObj) 
props <- props[ , c("feature_id", repr.Samps)]
colnames(props)[2:ncol(props)] <- cells

dtest.hcl <- dmTest(drimObj, coef = "CTSOXP")#contrast = c(0, 1))
saveRDS(dtest.hcl, file.path(opDir, paste(OP, "dobj.Min20.MF5.test", sep = ".")))

res.hcl <- results(dtest.hcl)
res.hcl$adj_pvalue[is.na(res.hcl$adj_pvalue)] <- 1
res.hcl.ft <- results(dtest.hcl, level = "feature")

######## Perform Stagewise analysis #########
resname <- "SoxNP.final.Min20.MF5.anu.out"
restest <- dtest.hcl
resdf <- res.hcl
resdf.ft <- res.hcl.ft

pConfirmation <- matrix(resdf.ft$pvalue, ncol = 1)
row.names(pConfirmation) <- gsub('\\.', 'p', resdf.ft$feature_id)
pScreen <- resdf$adj_pvalue
names(pScreen) <- gsub('\\.', 'p', resdf$gene_id)
  
tx2gene <- resdf.ft[ , c("feature_id", "gene_id")]
tx2gene$feature_id <- gsub('\\.', 'p', tx2gene$feature_id)  
tx2gene$gene_id <- gsub('\\.', 'p', tx2gene$gene_id)
  
sRObj.clcm <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted=TRUE, tx2gene=tx2gene)
sRObj.clcm <- stageWiseAdjustment(object = sRObj.clcm, method="dtu", alpha=0.05)

padj <- getAdjustedPValues(sRObj.clcm, order=TRUE, onlySignificantGenes=FALSE)
padj$geneID <- gsub("p", "\\.", padj$geneID)
padj$txID <- gsub("p", "\\.", padj$txID)
colnames(padj)[1:2] <- c("EID", "feature_id") 

final <- right_join(e2g, padj, by = "EID")
final[ ,c("ISOS", "SPLCHAIN", "OFLANK.5P", "OFLANK.3P")] <- NULL
final <- left_join(final, props, by = "feature_id")
final <- left_join(final, txcounts, by = "feature_id")
final <- left_join(final, nrcounts, by = "EID")

write.table(final, file.path(opDir, resname), quote = F, sep = "\t", row.names = F)







