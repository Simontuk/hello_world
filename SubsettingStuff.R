########
## Using new dataset, reading in with fread
########
# TODO: ADAPT TO command args!
library(minfi)
library(data.table)
library(MethylAid)



dataset <- "GBM"
mRNAseq <- paste(dataset,'.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt',sep="")
mRNAdata <- file.path("/Volumes/DEBINST/Data_GBM/gdac.broadinstitute.org_GBM.mRNAseq_Preprocess.Level_3.2015020400.0.0/GBM.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt")
data450k <- file.path("/Volumes/DEBINST/Data_GBM/gdac.broadinstitute.org_GBM.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0/GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")

# Source the Function file
source("~/Bachelor-Arbeit/Bachelor-Arbeit/scripts/00_GeneralFunctions.R")
source('~/hello_world/00_NewFunctions.R')

X.gbm <- fread(mRNAdata)
M.gbm <- readTCGA(data450k)

# Changing Column Names / Sample Names of Methylation Data Set
M.id <- colnames(M.gbm)
M.id <- sapply(lapply(strsplit(M.id,"\\-"),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})
colnames(M.gbm) <- M.id

# Changing Column Names / Sample Names of Expression Data Set
X.gbm$gene <- sapply(strsplit(X.gbm$gene,"\\|"),'[[',2)
setkey(X.gbm,gene)

#X.gbm <- X.gbm[,!"gene",with = F]
X.id <- colnames(X.gbm[,.SD,.SDcols=-key(X.gbm)])
X.id <- sapply(lapply(strsplit(X.id,"\\-"),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})

# Keep in Mind: data.table object
setnames(X.gbm,colnames(X.gbm[,.SD,.SDcols=-key(X.gbm)]),new = X.id)

# Gene-id extract

i.na <- which(apply(X.gbm[,.SD,.SDcols=-key(X.gbm)],1,function(x) {sum(!is.na(x))==0}))
X.entity <- X.gbm[-i.na,]
#dim(X.entity)


i.na <- which(apply(getBeta(M.gbm),1,function(x) {sum(!is.na(x))==0}))
M.entity <- M.gbm[-i.na,]
#dim(M.entity)


common <- intersect(X.id,M.id)
i.m <- match(common,M.id)
i.x <- match(common,X.id)

X.entity <- X.entity[,c(1,i.x+1),with=F]
M.entity <- M.entity[,i.m]

i.cg <- grep("cg",rownames(M.entity))
M.entity <- M.entity[i.cg,]


message('Computing median and sd ....')
STAT.entity <- data.frame(sd=apply(getBeta(M.entity),1,sd,na.rm=T),
                          iqr=apply(getBeta(M.entity),1,IQR,na.rm=T),
                          med=apply(getBeta(M.entity),1,median,na.rm=T),
                          mean=apply(getBeta(M.entity),1,mean,na.rm=T))

STAT.entity <- STAT.entity[order(STAT.entity$iqr,decreasing=T),]

####
# Choose the 100k most variable
####
#selected.cg <- rownames(STAT.entity)[1:100000]
#M.entity <- M.entity[selected.cg,]
####


####
#Saving the Data
####
# save('X.entity',
#      'M.entity',
#      'STAT.entity',
#      file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")),
#      compress='bzip2')
####
# rm()
# gc()
####
library(minfi)
library(data.table)
library(MethylAid)
library(parallel)
library(limma)
library(FDb.InfiniumMethylation.hg19)         ## Load Illumina 450k array
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

load("~/Bachelor-Arbeit/Bachelor-Arbeit/roadmap.seg.enh.RData")
roadmap.seg.enh <- GRangesList(roadmap.seg.enh)

### Data Processing



gid <- X.entity[,gene]
X.entity <- normalizeQuantiles(X.entity[,.SD,.SDcols=-key(X.entity)])

ALL <- colnames(X.entity)

tid <- sapply(strsplit(ALL,"\\."),'[',4)

TUMOR <- ALL[tid=='01']
CTRL <- ALL[tid=='11']


tx.red <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene"))
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes <- unlist(tx.red)
genes$gene_id <- names(genes)

# Different Subsets:
# Promoter:
prom.2kb <- promoters(genes,upstream=2000,downstream=500)
#Promoter from transcripts
prom.tx.2kb <- promoters(tx,upstream=2000,downstream=0)
#Genebody as defined by start-region (why?)
genebody.2kb <- promoters(genes,upstream=0,downstream=2000)

##Define TSS
tss <- promoters(genes,upstream=1,downstream=0)

## Promoter (2000 upstream, 500 downstream)
## SELECT MOST VARIABLE CpGs IN DIFFERENT GENOMIC COMPARTMENTS
selected.cg <- rownames(M.entity)

## promoter CpGs
ov <- findOverlaps(granges(M.entity[selected.cg]),prom.2kb)
i.prom.cg <- unique(queryHits(ov))
selected.prom <- selected.cg[i.prom.cg]

## gene-body CpGs
ov <- findOverlaps(granges(M.entity[selected.cg]),genes)
i.gb.cg <- unique(queryHits(ov))
selected.gb <- selected.cg[setdiff(i.gb.cg,i.prom.cg)]

## distal CpGs
selected.dist <- selected.cg[-c(i.prom.cg,i.gb.cg)]



cor.prom <- corMX1(selected.prom,w=0,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL')
cor.prom.rnd.gene <- corMX1(selected.prom,w=0,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL',randproc='gene')
cor.prom.rnd.sh <- corMX1(selected.prom,w=0,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL',randproc='samples')


cor.gb <- corMX1(selected.gb,w=0,X.entity,M.entity,gid,prom=genebody.2kb,subset='ALL')
cor.gb.rnd.gene <- corMX1(selected.gb,w=0,X.entity,M.entity,gid,prom=genebody.2kb,subset='ALL',randproc='gene')
cor.gb.rnd.sh <- corMX1(selected.gb,w=0,X.entity,M.entity,gid,prom=genebody.2kb,subset='ALL',randproc='samples')

##Selecting GB Enhancers
#CG.GENE <- cor.gb$CG
CG.gr <- granges(M.entity)

CGenh.gr <- subsetByOverlaps(CG.gr,roadmap.seg.enh)

cor.gb.enh <- cor.gb[cor.gb$CG %in% names(CGenh.gr),]
cor.gb.enh.rnd.gene <- cor.gb.rnd.gene[cor.gb.rnd.gene$CG %in% names(CGenh.gr),]
cor.gb.enh.rnd.sh <- cor.gb.rnd.sh[cor.gb.rnd.sh$CG %in% names(CGenh.gr),]
##

w <- 50000
cor.dist <- corMX1(selected.dist,w=w,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL')
cor.dist.rnd.gene <- corMX1(selected.dist,w=w,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL',randproc='gene')
cor.dist.rnd.sh <- corMX1(selected.dist,w=w,X.entity,M.entity,gid,prom=prom.2kb,subset='ALL',randproc='samples')

##Selecting Distal Enhancers
CG.GENE <- cor.dist$CG
CG.gr <- hm450[CG.GENE]

CGenh.gr <- subsetByOverlaps(CG.gr,roadmap.seg.enh)

cor.dist.enh <- cor.dist[cor.dist$CG %in% names(CGenh.gr),]
cor.dist.enh.rnd.gene <- cor.dist.rnd.gene[cor.dist.rnd.gene$CG %in% names(CGenh.gr),]
cor.dist.enh.rnd.sh <- cor.dist.rnd.sh[cor.dist.rnd.sh$CG %in% names(CGenh.gr),]
##



###

###