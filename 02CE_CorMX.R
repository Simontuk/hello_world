library(parallel)
library(limma)
library(FDb.InfiniumMethylation.hg19)         ## Load Illumina 450k array
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)


hm450 <- get450k()
load("/icgc/dkfzlsdf/analysis/B080/steiger/data/roadmap/roadmap.seg.enh.RData")
roadmap.seg.enh <- GRangesList(roadmap.seg.enh)


args <- commandArgs(TRUE)


dataset1 <- args[1]
dataset2 <- args[2]

dataset <- paste(dataset1,dataset2,sep="_")
#dataset1 <- "UCS"
#dataset2 <- "ACC"
## adapt path to dataset
main.dir <- '/icgc/dkfzlsdf/analysis/B080/steiger'
out.dir <- file.path(main.dir,"analysis","Cross-Entity",dataset)
data.dir <- file.path(out.dir,"RData")



r.dir <- file.path(main.dir,"scripts")
fig.dir <- file.path(out.dir,"figures")

source(file.path(r.dir,'00_GeneralFunctions.R'))

#Reading RDS (parallelized Dataset for every single Variable)
load(file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")))

X.entity <- normalizeQuantiles(X.entity)

ALL <- colnames(X.entity)

tid <- sapply(strsplit(ALL,"\\."),'[',5)

pid <- sapply(strsplit(ALL,"\\."),'[',1)

TUMOR <- ALL[tid=='01']
CTRL <- ALL[tid=='11']
dat1 <- ALL[pid==dataset1]
dat2 <- ALL[pid==dataset2]

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
     file=file.path(out.dir,paste('RData/',dataset,'.RData',sep="")))
