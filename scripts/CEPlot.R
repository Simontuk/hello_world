library(TeachingDemos)

args <- commandArgs(TRUE)

dataset1 <- args[1]
dataset2 <- args[2]

dataset <- paste(dataset1,dataset2,sep="_")


## adapt path to dataset
main.dir <- '/icgc/dkfzlsdf/analysis/B080/steiger'
out.dir <- file.path(main.dir,"analysis","Cross-Entity",dataset)
data.dir <- file.path(out.dir,"RData")

meth.dir <- file.path(data.dir)
file.dir <- paste('gdac.broadinstitute.org_',dataset,'.mRNAseq_Preprocess.Level_3.2015020400.0.0',sep="")
exp.dir <- file.path(data.dir,file.dir)

r.dir <- file.path(main.dir,"scripts")
fig.dir <- file.path(main.dir,"analysis/figures/Cross-Entity")

load(file.path(out.dir,paste('RData/',dataset,'.RData',sep="")))
load(file=file.path(out.dir,paste('RData/',dataset,'.load.RData',sep="")))

x <- function(){plot(density(STAT.entity[,]$mean),main=" Mean CpG-Methylation",ylab="",yaxt="n",xaxt="n",xlab="");
  lines(density(STAT.entity[rownames(M.entity),]$mean),col=5);
  lines(density(STAT.entity[cor.prom$CG,]$mean),col="red")
}
png(file.path(fig.dir,paste(sep="",dataset,"_density.promoters.png")))
par(mar=c(2,1,2,1))
plot(density(cor.prom.rnd.gene$pearson,na.rm = T),main=paste(dataset,"Density Promoters"),xlim=c(-1,1),yaxt='n',ylab=" ")
lines(density(cor.prom$pearson,na.rm = T),col="red")
legend("topright",legend = c("Promoters","Random Promoters"),
       col =  c("red","black"),pch=16)
legend("topleft",legend = c(paste("n = ",sum(!is.na(cor.prom.rnd.gene$pearson)))),
       col = c("red"), pch = 16)
subplot(x(),
  grconvertX(c(0.7, 0.95), "npc"), grconvertY(c(0.5,0.7), "npc"))

dev.off()



x <- function(){plot(density(STAT.entity[,]$mean),main=" Mean CpG-Methylation",ylab="",yaxt="n",xaxt="n",xlab="");
  lines(density(STAT.entity[rownames(M.entity),]$mean),col=5);
  lines(density(STAT.entity[cor.dist$CG,]$mean),col="red");
  lines(density(STAT.entity[cor.dist.enh$CG,]$mean),col="green")
}
png(file.path(fig.dir,paste(sep="",dataset,"_density.distal.png")))
par(mar=c(2,1,2,1))
plot(density(cor.dist.rnd.gene$pearson,na.rm = T),main=paste(dataset,"Density Distal"),xlim=c(-1,1),yaxt='n',ylab=" ")
lines(density(cor.dist.enh$pearson,na.rm = T),col="green")
lines(density(cor.dist$pearson,na.rm = T),col="red")
lines(density(cor.dist.enh.rnd.gene$pearson,na.rm = T),col="grey")
legend("topright",legend = c("Distal","Distal Enhancer","Random Distal","Random Dist-Enhancer"),
       col =  c("red","green","black","grey"),pch=16)
legend("topleft",legend =
       c(paste("n = ",sum(!is.na(cor.dist.rnd.gene$pearson))),
               paste("n = ",sum(!is.na(cor.dist.enh.rnd.gene$pearson)))),
       col = c("red","green"), pch = 16
       )
subplot(x(),
  grconvertX(c(0.7, 0.95), "npc"), grconvertY(c(0.5,0.7), "npc"))

dev.off()

x <- function(){plot(density(STAT.entity[,]$mean),main="Mean CpG-Methylation",ylab="",yaxt="n",xaxt="n",xlab="");
  lines(density(STAT.entity[rownames(M.entity),]$mean),col=5);
  lines(density(STAT.entity[cor.gb$CG,]$mean),col="red");
  lines(density(STAT.entity[cor.gb.enh$CG,]$mean),col="green")
}
png(file.path(fig.dir,paste(sep="",dataset,"_density.genebody.png")))
par(mar=c(2,1,2,1))
plot(density(cor.gb.rnd.gene$spearman,na.rm = T),main=paste(dataset,"Density Genebody"),xlim=c(-1,1),yaxt='n',ylab=" ")
lines(density(cor.gb.enh$spearman,na.rm = T),col="green")
lines(density(cor.gb$spearman,na.rm = T),col="red")
lines(density(cor.gb.enh.rnd.gene$spearman,na.rm = T),col="grey")
legend("topright",legend = c("Genebody","Genebody Enhancer","Random Genebody","Random GB-Enhancer"),
       col =  c("red","green","black","grey"),pch=16)
legend("topleft",legend =
       c(paste("n = ",sum(!is.na(cor.gb.rnd.gene$spearman))),
             paste("n = ",sum(!is.na(cor.gb.enh.rnd.gene$spearman)))),
       col = c("red","green"), pch = 16
       )
subplot(x(),
  grconvertX(c(0.7, 0.95), "npc"), grconvertY(c(0.5,0.7), "npc"))

dev.off()





#png(file.path(fig.dir,paste(sep="",dataset,"_distribMeth.png")))
#plot(density(M.entity[,1],na.rm = T),xlim=c(0,1),ylim=c(0,5))
#for (i in 1:ncol(M.entity)){
#  lines(density(M.entity[,i],na.rm = T))
#}
#dev.off()

#png(file.path(fig.dir,paste(sep="",dataset,"_MeanMeth.png")))
#plot(density(STAT.entity$mean),main=paste(dataset," Mean CpG-Methylation"),ylab="",yaxt="n")
#dev.off()
