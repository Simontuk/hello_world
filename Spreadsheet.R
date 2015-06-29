par(mfrow=c(1,1))
plot(as.numeric(getBeta(M.entity[35690])))
plot(as.numeric(XX[16778]))

PlotMX(as.numeric(getBeta(M.entity[117428])),(as.numeric(X.entity[2486])))



M.test <- readTCGA("/Volumes/DEBINST/Data_UCS/gdac.broadinstitute.org_UCS.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0/UCS.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
GBMtar.test <- readTCGA("/Volumes/DEBINST/Data_GBM/gdac.broadinstitute.org_GBM.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0.tar.gz")
?readTCGA

?tar
(M.test)

cross.test <- cbind(M.test,GBM.test)

X<-model.matrix(~M.test$Status)
minfi::getIslandStatus(minfi::makeGenomicRatioSetFromMatrix(as.matrix(M.entity)))

ALL <- colnames(cross.test)

tid <- sapply(strsplit(ALL,"\\-"),'[',4)
tid <- substr(tid,1,2)

cid <- sapply(strsplit(ALL,"\\-"),'[',2)
cid


UCS <- cbind(as.matrix(colnames(M.test)),1)
GBM <- cbind(as.matrix(colnames(GBM.test)),2)


ALL <- rbind(UCS,GBM)

indexUCS <- ALL[,2]==1
indexGBM <- ALL[,2]==2

X <- model.matrix(~ALL[,2])

require(doParallel)
registerDoParallel(cores = 6)
bumps <- bumphunter(cross.test,X,cutoff=00)
bumpn <- bumps$fitted[!is.na(bumps$fitted)]

plot(density(bumpn))

a <- names(bumps$fitted[which(bumps$fitted>0.1),1])
b <- rownames(M.entity)

pData(cross.test)

minfi::densityPlot(getBeta(cross.test),sampGroups = ALL[,2])
plotCpg(getBeta(cross.test),"cg16535999",pheno=ALL[,2],type = "cat")



sum(a %in% b)

