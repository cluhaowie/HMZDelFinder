##'---------------------------------------------------
##'Map candidate exon back to bedOrdered file
##'
##'
##' @param filtercandidateCalls    object return from callCandidateExon()
##' @param candidateZscore         object return from callCandidateExon()
##' @param bedOrdered
##' 
##'---------------------------------------------------
prepareExons <- function(filtercandidateCalls,bedOrdered,candidateZscore,mc.cores=4){
  if(is.null(names(filtercandidateCalls))){print("please check the names of filtercandidateCalls");return(NULL)}
  library(pbmcapply)
  library(data.table)
  temCallfun <- function(i,filtercandidateCalls,bedOrdered,candidateZscore,n=names(filtercandidateCalls)){
    gc()
    callname <- n[i]
    idx <- as.numeric(unlist(filtercandidateCalls[i]))
    data.table(bedOrdered[idx,],V5=candidateZscore[idx,callname])
  }
  print("[******Preparing DEL and DUP calls******]")
  temp <- pbmclapply(seq_along(filtercandidateCalls),temCallfun,filtercandidateCalls,bedOrdered,candidateZscore,mc.cores = mc.cores,max.vector.size =6656)
  names(temp)<- names(filtercandidateCalls)
  candidateExons <- rbindlist(temp,idcol=TRUE)
  colnames(candidateExons)[1] <- "Sample"
  rm(temp);gc()
  return(candidateExons) 
}

##'---------------------------------------------------
##'Merge exons to make a consistent call
##'
##' @param candiateExonCalls   object return from prepareExons
##' @param bedOrdered          ordered bed file
##' @param maxGap              The maxGap can be tolerate by the algorithm
##'----------------------------------------------------
mergeCandidates <- function(candidateExonsCalls, bedOrdered, maxGap = 10,mc.cores=2){
  resTmplist <- split.data.frame(candidateExonsCalls,candidateExonsCalls$Sample)
  MergeFUN <- function(x){
    maxGap <- 10
    x$mark_num<-1; x$exon_num<-1
    x$key <-paste(x$V1,"_",x$V2, sep="")
    bedOrderedIdx<- match(x$key, paste(bedOrdered$V1,"_", bedOrdered$V2, sep=""))
    df <- data.frame()
    j <- 1
    v  <-  bedOrdered[bedOrderedIdx[1],]
    v$mark_num<- 1;v$exon_num <- 1
    v$start_idx<- bedOrderedIdx[1]
    v$Zscore <- round(x$V5[1],3)
    diffIdx <- diff(bedOrderedIdx)
    for (i in diffIdx){
      j <- j + 1
      if (i >= maxGap){
        df <- rbindlist(list(df, v))
        v <-  bedOrdered[bedOrderedIdx[j],]
        v$mark_num <- 1
        v$exon_num <- 1
        v$start_idx<- bedOrderedIdx[j]
        v$Zscore <- round(x$V5[j],3)
      }else{
        v$mark_num <- v$mark_num+1 
        v$V3 <-bedOrdered[bedOrderedIdx[j],"V3"] # replace stop
        gn <- bedOrdered[bedOrderedIdx[j],"V4"] # gene name
        Zscore <- round(x$V5[j],3)
        if (gdata::trim(gn) != "" && !(gn %in% strsplit(v$V4,",")[[1]])){ v$V4 <- paste(v$V4, gn, sep=",")}
        v$Zscore <- paste(v$Zscore, Zscore, sep=",")
        v$exon_num <- bedOrderedIdx[j] - v$start_idx + 1
        z <- bedOrdered[bedOrderedIdx[j],]
      }}
    df <- rbindlist(list(df, v))
    df$Sample<-x$Sample[1]
    x <- as.data.frame(df)
    data.table(x[,c("V1", "V2", "V3", "V4","Zscore","start_idx", "mark_num","exon_num","Sample")])
  }
  resTmp <- pbmclapply(resTmplist,MergeFUN,mc.cores = mc.cores)
  res <- rbindlist(resTmp)
  rm(resTmplist,resTmp);gc()
  res <- res[order(res$V1, res$V2),]	
  colnames(res) <- c("Chr", "Start","Stop","Genes","Zscore","Start_idx","Mark_num","Exon_num","FID")
  res$Length <- res$Stop - res$Start
  res[order(res$Length, decreasing=T), ]
  res
}
