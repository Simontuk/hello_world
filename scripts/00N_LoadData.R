library(minfi)
library(data.table)
library(MethylAid)
library(parallel)
library(data.table)

args <- commandArgs(TRUE)


dataset <- args[1]
dataset <- "ACC"

## adapt path to dataset
main.dir <- '/icgc/dkfzlsdf/analysis/B080/steiger'
data.dir <- file.path(main.dir,'data',dataset)
out.dir <- file.path(data.dir,'output')
meth.dir <- file.path(data.dir,paste('gdac.broadinstitute.org_',dataset,'.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0',sep=""))
RNAseq.dir <- file.path(data.dir,paste('gdac.broadinstitute.org_',dataset,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep=""))
file.dir <- paste('gdac.broadinstitute.org_',dataset,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp.dir <- file.path(data.dir,file.dir)
r.dir <- file.path(main.dir,"scripts")
dir.create(out.dir)
dir.create(paste(out.dir,"/RData",sep=""))

mRNAseq <- paste(dataset,'.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt',sep="")
mRNAdata <- file.path(RNAseq.dir,mRNAseq)
data450k <- file.path(meth.dir,paste(dataset,'.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt',sep=""))

# Source the Function file
source("~/Bachelor-Arbeit/Bachelor-Arbeit/scripts/00_GeneralFunctions.R")
source('~/hello_world/00_NewFunctions.R')

X.entity <- fread(mRNAdata)
X.entity <- readTCGA(data450k)

# Changing Column Names / Sample Names of Methylation Data Set
M.id <- colnames(M.entity)
M.id <- sapply(lapply(strsplit(M.id,"\\-"),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})
colnames(M.entity) <- M.id

# Changing Column Names / Sample Names of Expression Data Set
X.entity$gene <- sapply(strsplit(X.entity$gene,"\\|"),'[[',2)
setkey(X.entity,gene)

#X.entity <- X.entity[,!"gene",with = F]
X.id <- colnames(X.entity[,.SD,.SDcols=-key(X.entity)])
X.id <- sapply(lapply(strsplit(X.id,"\\-"),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})

# Keep in Mind: data.table object
setnames(X.entity,colnames(X.entity[,.SD,.SDcols=-key(X.entity)]),new = X.id)

# Gene-id extract

i.na <- which(apply(X.entity[,.SD,.SDcols=-key(X.entity)],1,function(x) {sum(!is.na(x))==0}))
X.entity <- X.entity[-i.na,]
#dim(X.entity)


i.na <- which(apply(getBeta(M.entity),1,function(x) {sum(!is.na(x))==0}))
M.entity <- M.entity[-i.na,]
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


###
Saving the Data
###
save('X.entity',
     'M.entity',
     'STAT.entity',
     file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")),
     compress='bzip2')
###
rm()
gc()
###