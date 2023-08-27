getNrCounts <- function(fcounts, pdata, eidcol, grpsize, domin){
  nrdf1 <- data.frame(EID = unique(fcounts[, eidcol]), stringsAsFactors = F)
  nrdf2 <- data.frame(EID = unique(fcounts[, eidcol]), stringsAsFactors = F)
  for(ct in unique(pdata$CellType)){
    cts <- pdata$Sample[pdata$CellType == ct]
    for(ei in 1:nrow(nrdf1)){
      eid <- nrdf1$EID[ei]
      nrs <- colSums(fcounts[fcounts[, eidcol] == eid, cts], na.rm = T)
      nrs <- nrs[order(-nrs)]
      ##print(nrs)
      nrdf1[ei, ct] <- median(nrs)
      if(domin == T){nrdf1[ei, ct] <- min(nrs[1:grpsize])}
      nrdf2[ei, ct] <- paste(nrs, collapse = ", ")
    }
  }
  return(list(nrdf1, nrdf2))
}

getNrCountsSampW <- function(fcounts, pdata, eidcol){
  nrdf <- data.frame(EID = unique(fcounts[, eidcol]), stringsAsFactors = F)
  for(sa in unique(pdata$Sample)){
    for(ei in 1:nrow(nrdf)){
      eid <- nrdf$EID[ei]
      nrs <- sum(fcounts[fcounts[, eidcol] == eid, sa], na.rm = T)
      ##print(nrs)
      nrdf[ei, sa] <- nrs ## paste(nrs, collapse = ", ")
    }
  }
  return(nrdf)
}

