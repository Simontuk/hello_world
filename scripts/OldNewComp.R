library(parallel)
library(limma)
library(FDb.InfiniumMethylation.hg19)         ## Load Illumina 450k array
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)


hm450 <- get450k()
load("/icgc/dkfzlsdf/analysis/B080/steiger/data/roadmap/roadmap.seg.enh.RData")
roadmap.seg.enh <- GRangesList(roadmap.seg.enh)


args <- commandArgs(TRUE)


dataset <- args[1]

## adapt path to dataset
main.dir <- '/icgc/dkfzlsdf/analysis/B080/steiger'
data.dir <- file.path(main.dir,'data',dataset)
out.dir <- file.path(data.dir,'output')
meth.dir <- file.path(data.dir)
file.dir <- paste('gdac.broadinstitute.org_',dataset,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp.dir <- file.path(data.dir,file.dir)
r.dir <- file.path(main.dir,"scripts")
fig.dir <- file.path(out.dir,"figures")

source(file.path(r.dir,'00_GeneralFunctions.R'))
asd.dir <- "/icgc/dkfzlsdf/analysis/B080/steiger/analysis/output/RData/"
#Reading RDS (parallelized Dataset for every single Variable)
X.entity.old <- readRDS(paste(asd.dir,"X.",dataset,".rds",sep=""))
M.entity.old <- readRDS(paste(asd.dir,"M.",dataset,".rds",sep=""))
STAT.entity.old <- readRDS(paste(asd.dir,"STAT.",dataset,".rds",sep=""))

load(file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")))

X.entity <- normalizeQuantiles(X.entity)


ALL <- colnames(X.entity)

tid <- sapply(strsplit(ALL,"\\."),'[',4)

TUMOR <- ALL[tid=='01']
CTRL <- ALL[tid=='11']

## build correspondance table txid <-> geneid
#gid <- sapply(strsplit(rownames(X.entity),"\\|"),'[[',2)
#gid <- rownames(X.entity)
#symb <- sapply(strsplit(rownames(X.entity),"\\|"),'[[',1)
#rownames(X.entity) <- gid

# Define gene promoters from transcripts
tx.red <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene"))
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes <- unlist(tx.red)
genes$gene_id <- names(genes)

prom.2kb <- promoters(genes,upstream=2000,downstream=0)
prom.tx.2kb <- promoters(tx,upstream=2000,downstream=0)
genebody.2kb <- promoters(genes,upstream=0,downstream=2000)

tss <- promoters(genes,upstream=1,downstream=0)

## SELECT MOST VARIABLE CpGs IN DIFFERENT GENOMIC COMPARTMENTS
selected.cg <- rownames(M.entity)

## promoter CpGs
ov <- findOverlaps(hm450[selected.cg],prom.2kb)
i.prom.cg <- unique(queryHits(ov))
selected.prom <- selected.cg[i.prom.cg]

## gene-body CpGs
ov <- findOverlaps(hm450[selected.cg],genebody.2kb)
i.gb.cg <- unique(queryHits(ov))
selected.gb <- selected.cg[i.gb.cg]

## distal CpGs
selected.dist <- selected.cg[-c(i.prom.cg,i.gb.cg)]

cor.prom <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL')
cor.prom.rnd.gene <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='gene')
cor.prom.rnd.sh <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='samples')


cor.gb <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genebody.2kb,subset='ALL')
cor.gb.rnd.gene <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genebody.2kb,subset='ALL',randproc='gene')
cor.gb.rnd.sh <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genebody.2kb,subset='ALL',randproc='samples')

##Selecting GB Enhancers
CG.GENE <- cor.gb$CG
CG.gr <- hm450[CG.GENE]

CGenh.gr <- subsetByOverlaps(CG.gr,roadmap.seg.enh)

cor.gb.enh <- cor.gb[cor.gb$CG %in% names(CGenh.gr),]
cor.gb.enh.rnd.gene <- cor.gb.rnd.gene[cor.gb.rnd.gene$CG %in% names(CGenh.gr),]
cor.gb.enh.rnd.sh <- cor.gb.rnd.sh[cor.gb.rnd.sh$CG %in% names(CGenh.gr),]
##

w <- 50000
cor.dist <- corMX(selected.dist,w=w,X.entity,M.entity,prom=prom.2kb,subset='ALL')
cor.dist.rnd.gene <- corMX(selected.dist,w=w,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='gene')
cor.dist.rnd.sh <- corMX(selected.dist,w=w,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='samples')

##Selecting Distal Enhancers
CG.GENE <- cor.dist$CG
CG.gr <- hm450[CG.GENE]

CGenh.gr <- subsetByOverlaps(CG.gr,roadmap.seg.enh)

cor.dist.enh <- cor.dist[cor.dist$CG %in% names(CGenh.gr),]
cor.dist.enh.rnd.gene <- cor.dist.rnd.gene[cor.dist.rnd.gene$CG %in% names(CGenh.gr),]
cor.dist.enh.rnd.sh <- cor.dist.rnd.sh[cor.dist.rnd.sh$CG %in% names(CGenh.gr),]
##



which.cor <- 'spearman'

dir.create(fig.dir)
jpeg(file.path(fig.dir,paste(sep="",dataset,"_cumul.promoter.jpeg")))
PlotCumDist(list(Observed=cor.prom,RandGene=cor.prom.rnd.gene,ShuffleSamples=cor.prom.rnd.sh),
            main=paste(sep="","Promoter methylation ",dataset))
dev.off()

jpeg(file.path(fig.dir,paste(sep="",dataset,"_cumul.genebody.jpeg")))
PlotCumDist(list(Observed=cor.gb,RandGene=cor.gb.rnd.gene,ShuffleSamples=cor.gb.rnd.sh),
            main=paste(sep="","Gene body methylation ",dataset))
dev.off()

jpeg(file.path(fig.dir,paste(sep="",dataset,"_cumul.distal.jpeg")))
PlotCumDist(list(Observed=cor.dist,RandGene=cor.dist.rnd.gene,ShuffleSamples=cor.dist.rnd.sh),
            main=paste(sep="","Distal methylation ",dataset))
dev.off()

jpeg(file.path(fig.dir,paste(sep="",dataset,"_cumul.GB_enhancer.jpeg")))
PlotCumDist(list(Observed=cor.gb.enh,RandGene=cor.gb.enh.rnd.gene,ShuffleSamples=cor.gb.enh.rnd.sh),
            main=paste(sep="","Enhancer methylation ",dataset))
dev.off()

jpeg(file.path(fig.dir,paste(sep="",dataset,"_cumul.DIST_enhancer.jpeg")))
PlotCumDist(list(Observed=cor.dist.enh,RandGene=cor.dist.enh.rnd.gene,ShuffleSamples=cor.dist.enh.rnd.sh),
            main=paste(sep="","Enhancer methylation ",dataset))
dev.off()


##GB and GB enhancers
jpeg(file.path(fig.dir,paste(sep="",dataset,"_hist.gb_enh.jpeg")))
hist(cor.gb.enh$spearman,prob=T,col=rgb(0,1,0,0.5),main=paste(sep="","Genebody vs Genebody+Enhancer ",dataset))
hist(cor.gb$spearman,prob=T,col=rgb(1,0,0,0.5),add=T)
hist(cor.gb.rnd.gene$spearman,prob=T,col=rgb(0,0,1,0.5),add=T)
dev.off()

##Dist and Dist enhancers
jpeg(file.path(fig.dir,paste(dataset,"_hist.dist_enh.jpeg")))
hist(cor.dist.enh$spearman,prob=T,col=rgb(0,1,0,0.5),main=paste("Distal vs Distal+Enhancer ",dataset))
hist(cor.dist$spearman,prob=T,col=rgb(1,0,0,0.5),add=T)
hist(cor.dist.rnd.gene$spearman,prob=T,col=rgb(0,0,1,0.5),add=T)
dev.off()

### some striking examples of correlation
#
## strong negative correlation
##PlotMX(M.entity[6257,],X.entity[6942,],subtypes=list(tumor=TUMOR,control=CTRL),main="GIPC2 / cg19766489")
#PlotMX(M.entity[27422,],X.entity[249,],main="GIPC2 / cg19766489")
#
#PlotMX(M.entity[16044,],X.entity[1133,],main="GIPC2 / cg19766489")
#
#PlotMXrname(M.acc,X.acc,cor.prom,rowname=2362)
#
## strong positive correlation
##PlotMX(M.entity[24541,],X.entity[11510,],subtypes=list(tumor=TUMOR,control=CTRL),main="NDN / cg25061289")
#PlotMX(M.entity[24541,],X.entity[11510,],main="NDN / cg25061289")


#Build track

file <- file.path(out.dir,'high.neg.cor.dist.bed')
track.neg <- BuildTrack(cor.dist,cormin=-1,cormax=-0.6,tx=genes)
write.table("track name='NegCor' description='Negative correlations'",file=file,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
write.table(track.neg,file=file,append=TRUE,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

file <- file.path(out.dir,'high.pos.cor.dist.bed')
track.pos <- BuildTrack(cor.dist,cormin=0.6,cormax=1,tx=genes)
write.table("track name='PosCor' description='Positive correlations'",
            file=file,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
write.table(track.pos,file=file,append=TRUE,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

file <- file.path(out.dir,'high.neg.cor.prom.bed')
track.neg <- BuildTrack(cor.prom,cormin=-1,cormax=-0.6,tx=genes)
write.table("track name='NegCorProm' description='Negative correlations Prom'",file=file,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
write.table(track.neg,file=file,append=TRUE,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

file <- file.path(out.dir,'high.pos.cor.prom.bed')
track.pos <- BuildTrack(cor.prom,cormin=0.6,cormax=1,tx=genes)
write.table("track name='PosCorProm' description='Positive correlations Prom'",
            file=file,quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
write.table(track.pos,file=file,append=TRUE,
            quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)


#

save('cor.dist',
     'cor.dist.rnd.gene',
     'cor.dist.rnd.sh',
     'cor.prom',
     'cor.prom.rnd.gene',
     'cor.prom.rnd.sh',
     'cor.gb',
     'cor.gb.rnd.gene',
     'cor.gb.rnd.sh',
     'cor.gb.enh',
     'cor.gb.enh.rnd.gene',
     'cor.gb.enh.rnd.sh',
     'cor.dist.enh',
     'cor.dist.enh.rnd.gene',
     'cor.dist.enh.rnd.sh',
     file=file.path(out.dir,paste('RData/',dataset,'.RData')))
