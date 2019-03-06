##'------------------------------------------
##' Plot RPKM distribution on specific gene region and sample
##' 
##' @param sample       internal BAB number(character)
##' @param medelTrack     data.table that store the BAB number FID, and metaproject
##' @param gene           gene name
##' @param rpkmOrdered    object from 
##' @param bedOrdered
##' @param outputDir 
##' @example sample<-"BAB2992";gene<-"RAI1";plotDistribution(sample,medelTrack,gene,rpkmDtOrdered,bedOrdered,outputDir)
##'
##'--------------------------------------------
plotDistribution<-function(sample,medelTrack,gene,rpkmDtOrdered,bedOrdered,outputDir){
  if(!length(medelTrack[BAB==sample,FIDs])){print("BAB number don't have a FID.");return(NULL)}
  if(!(medelTrack[BAB==sample,FIDs] %in% rownames(rpkmDtOrdered))){print("FID don't have a rpkm record.");return(NULL)}
  if(!length(which(bedOrdered$V4==gene))){print("Gene is not in the bed file.");return(NULL)}
  sampleFID<-medelTrack[BAB==sample,FIDs]
  replaceInf <- function(x) {
    dm <- data.matrix(x)
    dm[!is.finite(dm)] <- 0
    res <- data.table(dm)
    res
  }
  window<-2
  trBlack <- rgb(0,0,0,alpha=0.3) 
  
  idx_lth<-which(bedOrdered$V4==gene)
  idx <- max(1, idx_lth[1]-window):min(idx_lth[length(idx_lth)]+window, ncol(rpkmDtOrdered))
  
  ll<-replaceInf(log(rpkmDtOrdered[,idx,with=F ] + 1, 10))
  rownames(ll)<-rownames(rpkmDtOrdered)
  png(paste0(outputDir,sample,"_",gene,".png",sep=""), width=1000,height=1000, pointsize=25)
  maxV <- max(t(ll))
  vspace <- 0.2* maxV
  par(  oma=c(2,2,3,2))
  alpha=0.5
  plot(c(min(idx), max(idx) ), c(-vspace, 1.05*maxV),type="n",xlim=c(min(idx), max(idx) ),xlab="Probe (exon) number",ylab="log(RPKM + 1)")
  matplot(matrix(rep(min(idx):max(idx), nrow(ll)), ncol=nrow(ll)), t(ll), type="l",lty=1, col=trBlack, add=T)
  lines(min(idx):max(idx),ll[which(rownames(rpkmDtOrdered)== sampleFID),],col="red",lwd=3)
  matplot(t(matrix(rep(idx,2),ncol=2)), t(cbind(rep(0, length(idx)),as.numeric(ll[which(rownames(ll)== sampleFID),]))),add=T,type="l", col="red",lwd=3, lty=1)
  abline(h=log(1+0.65, 10), lwd=3, lty=2, col="blue")
  abline(h=log(1+2*0.65, 10), lwd=3, lty=2, col="blue")
  
  ## plotting genes and exons
  geneIdx <- which(gene == bedOrdered$V4[idx])
  voffset <- -vspace * 0.15
  rect(idx[geneIdx]-0.2,rep(voffset + 0.03*vspace, length(idx[geneIdx])), idx[geneIdx] + 0.2, rep(voffset - 0.03*vspace, length(idx[geneIdx])), col="darkblue", border="darkblue")
  lines(c(min(idx[geneIdx]),max(idx[geneIdx])), c(voffset, voffset), col="darkblue", lwd=3)
  text(mean(idx[geneIdx]), voffset-0.16*vspace, gene, font=3, cex=0.7)
  dev.off()
}

##'------------------------------------------
##' Plot RPKM distribution on specific gene region and sample
##' 
##' @param sample          internal BAB number(character)
##' @param medelTrack      data.table that store the BAB number FID, and metaproject
##' @param gene            gene name
##' @param rpkmOrdered     object from 
##' @param bedOrdered
##' @param outputDir
##' @param distMat         distance matrix that match rowname of rpkmOrdered; 
##'				for example person correlation coefficiency matrix(log transformed is preferred) ; 
##'				 perClusterMat <- cor(t(rpkmDtOrdered[,sample(1:196907,10000)]));
##'				 rownames(perClusterMat)<-rownames(rpkmOrdered)
##' @example sample<-"BAB2992";gene<-"RAI1";plotDistribution(sample,medelTrack,gene,rpkmDtOrdered,bedOrdered,outputDir)
##'
##'--------------------------------------------
plotDistribution_TOP100 <- function(sample,medelTrack,gene,rpkmDtOrdered,bedOrdered,outputDir,distMat){
  if(!length(medelTrack[BAB==sample,FIDs])){print("BAB number don't have a FID.");return(NULL)}
  if(!(medelTrack[BAB==sample,FIDs] %in% rownames(rpkmDtOrdered))){print("FID don't have a rpkm record.");return(NULL)}
  if(!length(which(bedOrdered$V4==gene))){print("Gene is not in the bed file.");return(NULL)}
  sampleFID<-medelTrack[BAB==sample,FIDs]
  replaceInf <- function(x) {
    dm <- data.matrix(x)
    dm[!is.finite(dm)] <- 0
    res <- data.table(dm)
    rownames(res) <- rownames(x)
    res
  }
  window<-2
  trBlack <- rgb(0,0,0,alpha=0.3) 
  
  idx_lth<-which(bedOrdered$V4==gene)
  idx <- max(1, idx_lth[1]-window):min(idx_lth[length(idx_lth)]+window, ncol(rpkmDtOrdered))
  sorted_inx<-sort(distMat[sampleFID,],decreasing=T,index.return=TRUE)[['ix']][1:100]
  
  ll<-replaceInf(log(rpkmDtOrdered[sorted_inx,idx,with=F ] + 1, 10))
  rownames(ll)<-rownames(rpkmDtOrdered)[sorted_inx]
  png(paste0(outputDir,sample,"_",gene,'_',"Top100",".png",sep=""), width=1000,height=1000, pointsize=25)
  maxV <- max(t(ll))
  vspace <- 0.2* maxV
  par(  oma=c(2,2,3,2))
  alpha=0.5
  plot(c(min(idx), max(idx) ), c(-vspace, 1.05*maxV),type="n",xlim=c(min(idx), max(idx) ),xlab="Probe (exon) number",ylab="log(RPKM + 1)")
  matplot(matrix(rep(min(idx):max(idx), nrow(ll)), ncol=nrow(ll)), t(ll), type="l",lty=1, col=trBlack, add=T)
  lines(min(idx):max(idx),ll[which(rownames(ll)== sampleFID),],col="red",lwd=3)
  matplot(t(matrix(rep(idx,2),ncol=2)), t(cbind(rep(0, length(idx)),as.numeric(ll[which(rownames(ll)== sampleFID),]))),add=T,type="l", col="red",lwd=3, lty=1)
  abline(h=log(1+0.65, 10), lwd=3, lty=2, col="blue")
  abline(h=log(1+2*0.65, 10), lwd=3, lty=2, col="blue")
  
  ## plotting genes and exons
  geneIdx <- which(gene == bedOrdered$V4[idx])
  voffset <- -vspace * 0.15
  rect(idx[geneIdx]-0.2,rep(voffset + 0.03*vspace, length(idx[geneIdx])), idx[geneIdx] + 0.2, rep(voffset - 0.03*vspace, length(idx[geneIdx])), col="darkblue", border="darkblue")
  lines(c(min(idx[geneIdx]),max(idx[geneIdx])), c(voffset, voffset), col="darkblue", lwd=3)
  text(mean(idx[geneIdx]), voffset-0.16*vspace, gene, font=3, cex=0.7)
  dev.off()
}	     

##'------------------------------------------------------------------
##'call Candidate dup/del exon based on Z-score
##'
##'
##' @param disMatOrdered   distance/similarity matrix(assuming there will be function to prepare for this matrix)
##' @param group(in subFun)           each column from disMatOrdered object;distMat distance matrix or similarity matrix(pearson correlation coefficiency)
##'disMatOrdered <-sapply(rownames(distMat),function(x){sort(distMat[x,],decreasing=T,index.return=TRUE)[['ix']][1:100]})
##' @param rpkmDtOrdered   log transformed rpkmDtOrdered, note: row names match
##' @param cutoff          Zscore cutoff, default= +-2
##' 
##' @example mclapply(distMatOrdered,callCandidateExon,rpkmDtOrdered=rpkmDtOrdered,cutoff=2,mc.cores = 8 )
##'------------------------------------------------------------------
callCandidateExon <- function(distMatOrdered,rpkmDtOrdered,cutoff=2, mc.cores = 4){
  callCandidateExon<-function(group,rpkmDtOrdered,cutoff){
    gc();
    group<-as.numeric(unlist(group))
    ll<-rpkmDtOrdered[group,]
    l <- as.numeric(ll[1,])
    rownames(ll)<-rownames(rpkmDtOrdered)[group]
    Vmedian <- colMedians(as.matrix(ll))
    Vsd <- colSds(as.matrix(ll))
    Zscore <- as.numeric((l-Vmedian)/Vsd)
    candidatesInx <- which(Zscore>cutoff|Zscore< -cutoff);
    Candidates <- list(Zscore=Zscore,calls= candidatesInx)
    rm(ll)
    gc();
    return(Candidates)
  }
  print("[******Calculate Zscore and get candidate exons**********]")
  candidateExon <- pbmclapply(distMatOrdered,callCandidateExon,rpkmDtOrdered,cutoff=cutoff,mc.cores=mc.cores,max.vector.size =6656)  ##43.54676 secs
  return(candidateExon)
}     
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
    data.table(bedOrdered[idx,],V5=candidateZscore[idx,callname,with=FALSE])
  }
  print("[******Preparing DEL and DUP calls******]")
  temp <- pbmclapply(seq_along(filtercandidateCalls),temCallfun,filtercandidateCalls,bedOrdered,candidateZscore,mc.cores = mc.cores,max.vector.size =6656)
  names(temp)<- names(filtercandidateCalls)
  candidateExons <- rbindlist(temp,idcol=TRUE)
  colnames(candidateExons)[1] <- "Sample"
  colnames(candidateExons)[6] <- "V5"
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

##------------------------------------------------------------------------------
##' Finds overlap between calls and AOH regions
##' 
##' 
##' @param candidatesMergedAnnotated		object returned by annotateCandidates
##' @param extAOH							object returned by prepareAOHData
##' @param minAOHsize							min AOH size 
##' @param minAOHsig						AOH signal threshold
##' @param mc.cores							number of cores (see pbmclapply)
##------------------------------------------------------------------------------


annotateAOH <- function(candidatesMergedAnnotated, extAOH, minAOHsize, minAOHsig, mc.cores){
  aohSize<-minAOHsize
  if (is.null(extAOH)){
    candidatesMergedAnnotated[,paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- TRUE
    return (candidatesMergedAnnotated)
  }
  tmpRes <- pbmclapply((1:length(unique(candidatesMergedAnnotated$FID))), function(i){ fid <- unique(candidatesMergedAnnotated$FID)[i]
  tmpCand <- candidatesMergedAnnotated[candidatesMergedAnnotated$FID == fid,]
  tmpCand [,paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- FALSE
  starts <- tmpCand$Start
  stops <- tmpCand$Stop
  stops [which(starts > stops)] <- starts[which(starts > stops)]
  finalCandGR <- GRanges(tmpCand$Chr, IRanges(starts, stops))
  selAOH_tmp <- extAOH[extAOH$Name == fid & extAOH$Length > aohSize & extAOH$seg.mean > minAOHsig,]
  selAOH_tmp_gr <- GRanges(selAOH_tmp$chrom, IRanges(selAOH_tmp$loc.start.ext, selAOH_tmp$loc.end.ext))
  mm_tmp <- as.matrix(findOverlaps(finalCandGR, selAOH_tmp_gr))
  tmpCand[unique(mm_tmp[,1]), paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- TRUE
  tmpCand		
  },mc.cores=mc.cores)
  res <- rbindlist(tmpRes)
  res
  
}
