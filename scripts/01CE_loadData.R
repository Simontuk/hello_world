#############################################################################
### Loads methylation and expression datasets and returns a
### data frame for each of these data types
#############################################################################
library(parallel)
library(data.table)

args <- commandArgs(TRUE)


dataset1 <- args[1]
dataset2 <- args[2]

dataset <- paste(dataset1,dataset2,sep="_")

## adapt path to dataset
main.dir <- '/icgc/dkfzlsdf/analysis/B080/steiger'
data1.dir <- file.path(main.dir,'data',dataset1)
out1.dir <- file.path(data1.dir,'output')
file1.dir <- paste('gdac.broadinstitute.org_',dataset1,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp1.dir <- file.path(data1.dir,file1.dir)

data2.dir <- file.path(main.dir,'data',dataset2)
out2.dir <- file.path(data2.dir,'output')
file2.dir <- paste('gdac.broadinstitute.org_',dataset2,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp2.dir <- file.path(data2.dir,file2.dir)

r.dir <- file.path(main.dir,"scripts")

#dir.create(out.dir)
#dir.create(paste(out.dir,"/RData",sep=""))

#source(file.path(r.dir,'00_GeneralFunctions.R'))


## A . READ METHYLATION DATA for dataset1
message(paste("Reading meth-data ",dataset1))
Mheader1 <- read.table(file.path(data1.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                     sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")

M5rows1 <- read.table(file.path(data1.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                       sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")


classes <- sapply(M5rows1,class)
M1 <- read.table(file.path(data1.dir,'methylation450_NArm.txt'),
                sep="\t",header=F,row.names=1,fill=TRUE,skip=2,
                #nrows=1000,
                dec=".",comment.char = "",stringsAsFactors = F,
                colClasses=classes)


message("Done")
M1.id <- colnames(Mheader1)
M1.id <- sapply(lapply(strsplit(M1.id,"\\."),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})
M1.id <- paste(dataset1,M1.id,sep=".")
colnames(M1) <- M1.id


## A . READ METHYLATION DATA for dataset2
message(paste("Reading meth-data ",dataset2))
Mheader2 <- read.table(file.path(data2.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                       sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")

M5rows2 <- read.table(file.path(data2.dir,'methylation450_NArm.txt'),stringsAsFactors = F,
                      sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")


classes <- sapply(M5rows2,class)
M2 <- read.table(file.path(data2.dir,'methylation450_NArm.txt'),
                 sep="\t",header=F,row.names=1,fill=TRUE,skip=2,
                 #nrows=1000,
                 dec=".",comment.char = "",stringsAsFactors = F,
                 colClasses=classes)


message("Done")
M2.id <- colnames(Mheader2)
M2.id <- sapply(lapply(strsplit(M2.id,"\\."),'[',1:4),function(x) {y <- paste0(x,collapse="."); gsub("[A-Z]$","",y)})
M2.id <- paste(dataset2,M2.id,sep=".")
colnames(M2) <- M2.id


## B. READ EXPRESSION DATA for dataset 1


message("Reading Expression-Data for dataset1")
mRNAseq <- paste(dataset1,'.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt',sep="")

Xheader1 <- read.table(file.path(exp1.dir,mRNAseq),stringsAsFactors = F,
                       sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")
X5rows1 <- read.table(file.path(exp1.dir,mRNAseq),stringsAsFactors = F,
                     sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")

classes <- sapply(X5rows1,class)

X1 <- read.table(file.path(exp1.dir,mRNAseq),
                  sep="\t",header=F,row.names=1,fill=TRUE,skip=1,
                  #nrows=1000,
                  dec=".",comment.char = "",stringsAsFactors = F,
                  colClasses=classes)

colnames(X1) <- paste(dataset1,colnames(Xheader1),sep=".")



## Reformat Geneid with strsplit
gid <- sapply(strsplit(rownames(X1),"\\|"),'[[',2)
rownames(X1) <- gid


## B. READ EXPRESSION DATA for dataset 2
message("Reading Expression-Data for dataset2")
mRNAseq <- paste(dataset2,'.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt',sep="")

Xheader2 <- read.table(file.path(exp2.dir,mRNAseq),stringsAsFactors = F,
                       sep="\t",header=T,row.names=1,fill=TRUE,nrow=1,dec = ".")
X5rows2 <- read.table(file.path(exp2.dir,mRNAseq),stringsAsFactors = F,
                      sep="\t",header=F,row.names=1,skip=2,fill=TRUE,nrow=5,dec = ".")

classes <- sapply(X5rows2,class)

X2 <- read.table(file.path(exp2.dir,mRNAseq),
                 sep="\t",header=F,row.names=1,fill=TRUE,skip=1,
                 #nrows=1000,
                 dec=".",comment.char = "",stringsAsFactors = F,
                 colClasses=classes)

colnames(X2) <- paste(dataset2,colnames(Xheader2),sep=".")


## Reformat Geneid with strsplit
gid <- sapply(strsplit(rownames(X2),"\\|"),'[[',2)
rownames(X2) <- gid

int <- intersect(rownames(X1),rownames(X2))
i.x1 <- match(int,rownames(X1))
i.x2 <- match(int,rownames(X2))

X1 <- X1[i.x1,]
X2 <- X2[i.x2,]

int <- intersect(rownames(M1),rownames(M2))
i.m1 <- match(int,rownames(M1))
i.m2 <- match(int,rownames(M2))

M1 <- M1[i.m1,]
M2 <- M2[i.m2,]


M <- cbind(M1,M2)
X <- cbind(X1,X2)

M.id <- colnames(M)
X.id <- colnames(X)

## SELECTING COMMON SAMPLES

common <- intersect(X.id,M.id)
i.m <- match(common,M.id)
i.x <- match(common,X.id)

#Removing NAs in X
i.na <- which(apply(X,1,function(x) {sum(!is.na(x))==0}))
X <- X[-i.na,]


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


out.dir <- file.path(main.dir,"analysis","Cross-Entity",dataset)
dir.create(out.dir)
dir.create(file.path(out.dir,"RData"))

save('X.entity',
    'M.entity',
    'STAT.entity',
    file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")),
    compress='bzip2')

rm()
gc()
