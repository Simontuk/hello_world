####
#CorMX1 für data.table and GenomicRatioSets
####

corMX1 <- function(cg,               # vector of CpG identifiers
                  w=0,               # in case you want to allow CpGs to be a little off the region of interest
                  X,                 # data.table with expression data
                  M,                 # GenomicRatioSet with methylation data
                  gid,               # Gene Names for the data.table in X  
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
  gid <- gid
  
  P <- prom[which(prom$gene_id %in% gid)]  
  hm450.sel <- hm450[which(names(hm450) %in% cg)]
  
  
  ## Old Calculation of Overlaps
  #ov <- findOverlaps(hm450.sel+w,P)
  #names(hm450.sel)
  #queryHits(ov)
  #CG <- names(hm450.sel)[queryHits(ov)]
  
  
  #hmw <- hm450.sel+w
  
  #
  ov <- findOverlaps(granges(M)+w,P)
  G <-  P$gene_id[subjectHits(ov)]
  
  #
  MM <- M[queryHits(ov)]
  rnmm <- sub("\\..*","",rownames(MM))
  MM.chr <- as.character(seqnames(MM))
  
  XX <- X[match(G,gid)]
  #gid <-  gid[match(G,gid)]
  rnxx <- sub("\\..*","",gid)
  XX.chr <- as.character(seqnames(P[rnxx]))
  
  #gidlist <- cbind(gid,rownames(XX))
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
  i.x <- match(rnxx,gid)
  
  message("Compute correlation")
  ## compute correlation
  CC <- do.call("rbind",mclapply(1:nrow(MM),function(i) {
    cp <- cor(as.numeric(getBeta(MM[i,])),as.numeric(XX[i]),use="everything")
    cs <- cor(as.numeric(getBeta(MM[i,])),as.numeric(XX[i]),use="everything",method="spearman")
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
  },mc.cores=8))
  
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



###
#Plot MX 
###
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
