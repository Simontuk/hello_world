dataset <- "THYM"
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
load(file=file.path(out.dir,paste('RData/',dataset,'.RData',sep=" ")))

png(file.path(fig.dir,paste(sep="",dataset,"_distribMeth.png")))
plot(density(M.entity[,1],na.rm = T),xlim=c(0,1),ylim=c(0,5))
for (i in 1:ncol(M.entity)){
  lines(density(M.entity[,i],na.rm = T))
}
dev.off()

png(file.path(fig.dir,paste(sep="",dataset,"_distribMethold.png")))
plot(density(as.numeric(M.entity.old[,1]),na.rm = T),ylim=c(0,9e-06),col=2)
for (i in 1:ncol(M.entity.old)){
  a <- i
  lines(density(as.numeric(M.entity.old[,i]),na.rm = T),col=2)
}
dev.off()