corMX1 <- function(cg,                # vector of CpG identifiers
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
  #names(hm450.sel)
  #queryHits(ov)
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
  },mc.cores=6))
  
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
  #cor.dist <- cc[order(cc$spearman),]
}
