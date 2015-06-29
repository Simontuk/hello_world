## load libraries
library(gplots)
#library(multicore)
library(org.Hs.eg.db)
library(FDb.InfiniumMethylation.hg19)         ## Load Illumina 450k array
library(TxDb.Hsapiens.UCSC.hg19.knownGene)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot methylation and expression profiles
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotMX <- function(M,X,main="",col=NULL,add=FALSE,xlim=NULL,
                   ylim=NULL,subtypes=NULL,loc='bottomleft') {
  
  require(RColorBrewer)
  if (is.null(col) && is.null(subtypes)) {col <- 'blue'}
  
  if (!is.null(subtypes)) {
    st <- as.factor(unlist(lapply(names(subtypes),function(x) {
      y <- rep(x,length(subtypes[[x]]));names(y) <- subtypes[[x]];return(y)})))
    color <- brewer.pal(max(3,length(subtypes)),"Set1")
    col <- color[st[names(M)]]
  }
  
  if (is.null(xlim)) {xlim=c(0,1)}
  x <- as.numeric(X)
  m <- as.numeric(M)
  i.inf <- union(which(is.infinite(x)),which(is.infinite(m)))
  if (length(i.inf)>0)
  {
    x <- x[-i.inf]
    m <- m[-i.inf]
  }
  xmax <- max(x)
  xmin <- min(x)
  ytextpos <- xmin+0.5*(xmax-xmin)
  
  if (add) {
    points(m,x,
           col=col,pch=19)
  }
  else {
    plot(m,x,xlim=xlim,ylim=ylim,
         col=col,pch=19,
         main=main,
         xlab="Methylation",ylab="Expression")
    a <- lm(xr~mr,data.frame(xr=x,mr=m))  
    #a <- lm(xr~mr,data.frame(xr=rank(x),mr=rank(m)))
    abline(a,lwd=2,col=col)
    #message(paste0("a=",a))
    text(0.1,ytextpos,labels=sprintf("R=%.2f",cor(m,x,method="spearman")))
    if (!is.null(subtypes)) {legend(loc,
                                    as.character(names(subtypes)),
                                    col=color[as.factor(names(subtypes))],
                                    pch=19)
                             mm <- sapply(subtypes,function(y) {mean(as.numeric(M[y]),na.rm=TRUE)})
                             xx <- sapply(subtypes,function(y) {mean(as.numeric(X[y]),na.rm=TRUE)})
                             #message(mm)
                             abline(v=mm,col=color[as.factor(names(subtypes))],lty=2)
                             abline(h=xx,col=color[as.factor(names(subtypes))],lty=2)
                             
    }
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot methylation and expression profiles only with rowname of Cor-Table
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotMXrname <- function(M,X,VV,rowname,col=NULL,add=FALSE,xlim=NULL,
                   ylim=NULL,subtypes=NULL,loc='bottomleft') {
  
  require(RColorBrewer)
  if (is.null(col) && is.null(subtypes)) {col <- 'blue'}
  
  if (!is.null(subtypes)) {
    st <- as.factor(unlist(lapply(names(subtypes),function(x) {
      y <- rep(x,length(subtypes[[x]]));names(y) <- subtypes[[x]];return(y)})))
    color <- brewer.pal(max(3,length(subtypes)),"Set1")
    col <- color[st[names(M)]]
  }
  
  xind <- VV[rowname,"index.x"]
  mind <- VV[rowname,"index.m"]
  main <- paste(VV[rowname,"CG"],VV[rowname,"symb"])
  
  if (is.null(xlim)) {xlim=c(0,1)}
  x <- as.numeric(X[xind,])
  m <- as.numeric(M[mind,])
  i.inf <- union(which(is.infinite(x)),which(is.infinite(m)))
  if (length(i.inf)>0)
  {
    x <- x[-i.inf]
    m <- m[-i.inf]
  }
  xmax <- max(x)
  xmin <- min(x)
  ytextpos <- xmin+0.5*(xmax-xmin)
  
  if (add) {
    points(m,x,
           col=col,pch=19)
  }
  else {
    plot(m,x,xlim=xlim,ylim=ylim,
         col=col,pch=19,
         main=main,
         xlab="Methylation",ylab="Expression")
    a <- lm(xr~mr,data.frame(xr=x,mr=m))  
    #a <- lm(xr~mr,data.frame(xr=rank(x),mr=rank(m)))
    abline(a,lwd=2,col=col)
    #message(paste0("a=",a))
    text(0.1,ytextpos,labels=sprintf("R=%.2f",cor(m,x,method="spearman")))
    text(0.1,ytextpos+0.5,labels=sprintf("p-val=%.8f",VV[rowname,"p.value.pears"]))
    if (!is.null(subtypes)) {legend(loc,
                                    as.character(names(subtypes)),
                                    col=color[as.factor(names(subtypes))],
                                    pch=19)
                             mm <- sapply(subtypes,function(y) {mean(as.numeric(M[y]),na.rm=TRUE)})
                             xx <- sapply(subtypes,function(y) {mean(as.numeric(X[y]),na.rm=TRUE)})
                             #message(mm)
                             abline(v=mm,col=color[as.factor(names(subtypes))],lty=2)
                             abline(h=xx,col=color[as.factor(names(subtypes))],lty=2)
                             
    }
  }
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Compute statistics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makeStat <- function(M,TYP=TYPES) {
  
  STAT <- sapply(TYP,function(t) {
    data.frame(stat.mean=apply(M[,get(t)],1,mean),
               stat.sd=apply(M[,get(t)],1,sd),
               stat.med=apply(M[,get(t)],1,median),
               stat.iqr=apply(M[,get(t)],1,IQR)
    )
  },simplify=FALSE)
  return(STAT)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### FDR related stuff
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fdr <- function(p) {
  length(pval.rand[pval.rand <= p])/length(pval.dmr[pval.dmr <= p])
}

plotFDR <- function(min=3,max=10) {
  p <- 10^(-seq(min,max,by=0.2))
  plot(-log10(p),sapply(p,fdr),xlab="-log10(pval)",ylab="fdr",type="l",main="FDR",col="red",lwd=2)
  abline(h=1e-3,lty=2)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   CORRELATION METHYLATION EXPRESSION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# for a given set of CpGs, width and expression dataset, compute the correlation
# between methylation profile and gene expression
# returns a data frame with pearson and spearman correlation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corMX <- function(cg,                # vector of CpG identifiers
                  w=0,               # in case you want to allow CpGs to be a little off the region of interest
                  X,                 # data.frame with expression data
                  M,                 # data.frame with methylation data
                  randproc='none',   # randomization procedure to get NULL distribution
                  prom,              # GRanges object with promoters
                  subset=NULL) {     # if you want to use a subset of the samples (i.e. columns)
  
  if (is.null(subset)) {ALL <- colnames(X);subset <- 'ALL'}
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  
  x <- org.Hs.egSYMBOL
  geneid2symb.df <- as.data.frame(x[mappedkeys(x)])
  geneid2symb <- as.list(geneid2symb.df[,2])
  names(geneid2symb) <- geneid2symb.df[,1]
  gid <- rownames(X)
  
  P <- prom[which(prom$gene_id %in% gid)]  
  hm450.sel <- hm450[which(names(hm450) %in% cg)]
  
  ov <- findOverlaps(hm450.sel+w,P)
  names(hm450.sel)
  queryHits(ov)
  CG <- names(hm450.sel)[queryHits(ov)]
  G <-  P$gene_id[subjectHits(ov)]
  
  MM <- M[match(CG,rownames(M)),get(subset)]
  rnmm <- sub("\\..*","",rownames(MM))
  MM.chr <- as.character(seqnames(hm450[rnmm]))
  XX <- X[match(G,gid),get(subset)]
  rnxx <- sub("\\..*","",rownames(XX))
  XX.chr <- as.character(seqnames(P[rnxx]))
  
  
  ## make sure to reorder the samples ...
  i <- match(colnames(XX),colnames(MM))
  MM <- MM[,i]
  
  
  ## randomization procedures
  if (randproc == 'gene') { 
    XX <- X[sample(1:nrow(X),nrow(XX),replace=TRUE),get(subset)]
  }
  else if (randproc == 'cpg') {
    M.dist <- M[which(rownames(M) %in% cg.dist),];
    MM <- M.dist[sample(1:nrow(M.dist),nrow(MM),replace=TRUE),get(subset)]
  }  
  else if (randproc == 'both') {
    M.dist <- M[which(rownames(M) %in% cg.dist),];
    MM <- M.dist[sample(1:nrow(M.dist),nrow(MM),replace=TRUE),get(subset)];
    XX <- X[sample(1:nrow(X),nrow(XX),replace=TRUE),get(subset)]    
  }
  
  else if (randproc == 'samples') {MM <- MM[,sample(get(subset))]}
  
  else if (randproc == 'diffchrom') {
    i<- sapply(XX.chr,function(cx) {sample(which(MM.chr != cx),1)})
    MM <- MM[i,] 
  }
  i.m <- match(rnmm,rownames(M))
  i.x <- match(rnxx,rownames(X))
  
  message("Compute correlation")
  ## compute correlation
  CC <- do.call("rbind",mclapply(1:nrow(MM),function(i) {
    cp <- cor(as.numeric(MM[i,]),as.numeric(XX[i,]),use="everything")
    cs <- cor(as.numeric(MM[i,]),as.numeric(XX[i,]),use="everything",method="spearman")
    #     if(sum(is.na(MM[i,]))>sum(is.na(XX[i,]))){
#       npr <- ncol(MM)-sum(is.na(MM[i,]))  
#     } else {
#       npr <- ncol(XX)-sum(is.na(XX[i,]))
#     }
     pval <- NA
#     if(sum(is.na(XX[i,]))>length(MM)-2 | sum(is.na(MM[i,]))>length(XX)-2){
#       pval <- NA
#     }
#     else{
#     
#       pval <- cor.test(as.numeric(MM[i,]),as.numeric(XX[i,]),method="spearman",exact=FALSE)$p.value
#       
#     }
    
    
    c(cp,cs)
  },mc.cores=12))
  
  message("Matching names")
  
  ii <- match(rnxx,names(geneid2symb))
  a <- geneid2symb[ii]
  i <- which(unlist(lapply(a,is.null))==TRUE)
  names(a)[i] <- rnxx[i]
  a[i] <- NA
  
  cc <- data.frame(CG=rnmm,
                   geneid=rnxx,
                   pearson=CC[,1],
                   spearman=CC[,2],
                   symb=unlist(a),
                   index.m=i.m,
                   index.x=i.x
                   
                   ) 
  
  rownames(cc) <- NULL
  return(cc[order(cc$spearman),])
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###   Build UCSC track
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BuildTrack <- function(corr,cormin=-1,cormax=-.5,col="255,0,0",tx) {
  CG.GENE <- corr[which(corr$spearman < cormax & corr$spearman > cormin),]
  CG.gr <- hm450[CG.GENE$CG]
  GENE.gr <- tx[match(CG.GENE$geneid,tx$gene_id)]
  TSS <- promoters(GENE.gr,upstream=1,downstream=0)
  do.call("rbind",
          lapply(1:length(CG.gr), function(i) {
            a <- CG.gr[i]
            b <- TSS[i]
            corr.pair <- CG.GENE[i,]$spearman
            start.int <- ifelse(start(a) < start(b),start(a),start(b))
            strand <- ifelse(start(a) < start(b),'+','-')
            stop.int <- ifelse(end(a) > end(b),end(a),end(b))
            c(as.character(seqnames(a)),
              start.int,
              stop.int,
              paste0(names(a),"_",as.character(CG.GENE[i,'symb']),"_",sprintf("%.2f",corr.pair)),
              0,
              #              abs(CG.GENE[i,]$pearson)*1000,
              strand,
              start.int,
              stop.int,
              col,
              2,
              "1,1",
              as.character(paste(0,stop.int-1-start.int,sep=","))
            )
          }))
}  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%  Plot the cumulative distribution of correaltion values
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotCumDist <- function(cor.list,which.cor='spearman',main=NULL,ylim=NULL) {
  
  require(RColorBrewer)
  n <- length(cor.list)
  col <- brewer.pal(n,'Set1')
  par(mfrow=c(2,1))
  
  cor <- cor.list[[1]]
  xx <- sort(cor[,which.cor],decreasing=TRUE)
  plot(xx,1:length(xx),type="s",log="y",col=col[1],
       main=main,
       lwd=2,
       xlim = c(-1,1),
       ylim=ylim,
       xlab=paste0("Correlation (",which.cor,")"),
       ylab="Cumulative distribution")
  
  lapply(2:length(cor.list),function(i) {
    cor <- cor.list[[i]]
    xx <- sort(cor[,which.cor],decreasing=TRUE)
    points(xx,1:length(xx),type="s",col=col[i],lwd=2)
  })
  legend("bottomleft",names(cor.list),col=col,lty=1,lwd=2)
  
  cor <- cor.list[[1]]  
  xx <- sort(cor[,which.cor],decreasing=TRUE)
  plot(xx,rev(1:length(xx)),type="s",log="y",col=col[1],
       main=main,
       lwd=2,
       xlim = c(-1,1),
       ylim = ylim,
       xlab=paste0("Correlation (",which.cor,")"),
       ylab="Cumulative distribution")
  
  lapply(2:length(cor.list),function(i) {
    cor <- cor.list[[i]]
    xx <- sort(cor[,which.cor],decreasing=TRUE)
    points(xx,rev(1:length(xx)),type="s",col=col[i],lwd=2)
  })
  legend("bottomright",names(cor.list),col=col,lty=1,lwd=2)               
}

PlotDist <- function(cor.list,which.cor='spearman',main="",ylim=NULL) {
  
  require(RColorBrewer)
  n <- length(cor.list)
  col <- brewer.pal(n,'Set1')
  
  cor <- cor.list[[1]]  
  xx <- sort(cor[,which.cor],decreasing=TRUE)
  plot(density(xx),col=col[1],
       main=main,
       lwd=2,
       ylim=ylim,
       xlab=paste0("Correlation (",which.cor,")"))
  
  lapply(2:length(cor.list),function(i) {
    cor <- cor.list[[i]]
    xx <- sort(cor[,which.cor],decreasing=TRUE)
    lines(density(xx),col=col[i],lwd=2)
  })
  legend("topleft",names(cor.list),col=col,lty=1,lwd=2)
  
}
