#############################################################################
### Loads methylation and expression datasets and returns a
### data frame for each of these data types
#############################################################################
library(parallel)
library(data.table)

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
dir.create(out.dir)
dir.create(paste(out.dir,"/RData",sep=""))

#source(file.path(r.dir,'00_GeneralFunctions.R'))


## A . READ METHYLATION DATA
message(paste("Reading meth-data ",dataset))
Mheader <- read.table(file.path(meth.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                     sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")

M5rows <- read.table(file.path(meth.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                       sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")


classes <- sapply(M5rows,class)

M <- read.table(file.path(meth.dir,'methylation450_NArm.txt'),
                sep="\t",header=F,row.names=1,fill=TRUE,skip=2,
                #nrows=1000,
                dec=".",comment.char = "",stringsAsFactors = F,
                colClasses=classes)


message("Done")
M.id <- colnames(Mheader)
M.id <- sapply(lapply(strsplit(M.id,"\\."),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})
colnames(M) <- M.id
## B. READ EXPRESSION DATA


message("Reading Expression-Data")
mRNAseq <- paste(dataset,'.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt',sep="")

Xheader <- read.table(file.path(exp.dir,mRNAseq),stringsAsFactors = F,
                       sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")
X5rows <- read.table(file.path(exp.dir,mRNAseq),stringsAsFactors = F,
                     sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")


classes <- sapply(X5rows,class)

X <- read.table(file.path(exp.dir,mRNAseq),
                  sep="\t",header=F,row.names=1,fill=TRUE,skip=1,
                  #nrows=1000,
                  dec=".",comment.char = "",stringsAsFactors = F,
                  colClasses=classes)

colnames(X) <- colnames(Xheader)


i.na <- which(apply(X,1,function(x) {sum(!is.na(x))==0}))
X <- X[-i.na,]
X.id <- colnames(X)

## build correspondance table txid <-> geneid
gid <- sapply(strsplit(rownames(X),"\\|"),'[[',2)
rownames(X) <- gid


## SELECTING COMMON SAMPLES

common <- intersect(X.id,M.id)
i.m <- match(common,M.id)
i.x <- match(common,X.id)

X.entity <- X[,i.x]
M.entity <- M[,i.m]

i.cg <- grep("cg",rownames(M.entity))
M.entity <- M.entity[i.cg,]

message('Computing median and sd ....')
STAT.entity <- data.frame(sd=apply(M.entity,1,sd,na.rm=T),
                     iqr=apply(M.entity,1,IQR,na.rm=T),
                     med=apply(M.entity,1,median,na.rm=T),
                     mean=apply(M.entity,1,mean,na.rm=T))


## rank CpGs by decreasing variability (ost variable on top)
STAT.entity <- STAT.entity[order(STAT.entity$iqr,decreasing=T),]
selected.cg <- rownames(STAT.entity)[1:100000]
M.entity <- M.entity[selected.cg,]


# saveRDS(X.entity,paste(out.dir,"/RData/","X.",dataset,".rds",sep=""),compress='bzip2')
# saveRDS(M.entity,paste(out.dir,"/RData/","M.",dataset,".rds",sep=""),compress='bzip2')
# saveRDS(STAT.M,paste(out.dir,"/RData/","STAT.",dataset,".rds",sep=""),compress='bzip2')


save('X.entity',
    'M.entity',
    'STAT.entity',
    file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")),
    compress='bzip2')

rm()
gc()
