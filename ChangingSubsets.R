library(parallel)
library(limma)
library(FDb.InfiniumMethylation.hg19)         ## Load Illumina 450k array
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)


hm450 <- get450k()
load("/icgc/dkfzlsdf/analysis/B080/steiger/data/roadmap/roadmap.seg.enh.RData")
roadmap.seg.enh <- GRangesList(roadmap.seg.enh)

dataset <- "UCS"
main.dir <- file.path('~/Bachelor-Arbeit/Bachelor-Arbeit')
data.dir <- file.path(main.dir,'data',dataset)
out.dir <- file.path(data.dir,'output')
meth.dir <- file.path(data.dir)
file.dir <- paste('gdac.broadinstitute.org_',dataset,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp.dir <- file.path(data.dir,file.dir)
r.dir <- file.path(main.dir,"scripts")
fig.dir <- file.path(out.dir,"figures")

source(file.path(r.dir,'00_GeneralFunctions.R'))

load(file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")))
load('~/Bachelor-Arbeit/Bachelor-Arbeit/UCS/output/RData/UCS.load.RData')


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
genebody.2kb <- promoters(genes,upstream=0,downstream=6000)



tss <- promoters(genes,upstream=1,downstream=0)

## SELECT MOST VARIABLE CpGs IN DIFFERENT GENOMIC COMPARTMENTS
selected.cg <- rownames(M.entity)

## promoter CpGs
ov <- findOverlaps(hm450[selected.cg],prom.2kb)
i.prom.cg <- unique(queryHits(ov))
selected.prom <- selected.cg[i.prom.cg]


ov <- findOverlaps(hm450[selected.cg],genes)
i.gb.cg <- unique(queryHits(ov))
selected.gb <- selected.cg[setdiff(i.gb.cg,i.prom.cg)]


## distal CpGs
selected.dist <- selected.cg[-c(i.prom.cg,i.gb.cg)]

##Testing
cor.gb <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genes,subset='ALL')
cor.gb.rnd.gene <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genes,subset='ALL',randproc='gene')
cor.gb.rnd.sh <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genes,subset='ALL',randproc='samples')
###

cor.prom <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL')
cor.prom.rnd.gene <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='gene')
cor.prom.rnd.sh <- corMX(selected.prom,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='samples')


cor.gb1new <- corMX(selected.gb,w=0,X.entity,M.entity,prom=genebody.2kb,subset='ALL')
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
