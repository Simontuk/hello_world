png(file.path(fig.dir,paste(sep="",dataset,"_density.genebody_newdef3.png")))
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
dev.off()


png(file.path(fig.dir,paste(sep="",dataset,"_density.distal_newset3.png")))
par(mar=c(2,1,2,1))
plot(density(cor.dist.rnd.gene$spearman,na.rm = T),main=paste(dataset,"Density Distal"),xlim=c(-1,1),yaxt='n',ylab=" ")
lines(density(cor.dist.enh$spearman,na.rm = T),col="green")
lines(density(cor.dist$spearman,na.rm = T),col="red")
lines(density(cor.dist.enh.rnd.gene$spearman,na.rm = T),col="grey")
legend("topright",legend = c("Genebody","Genebody Enhancer","Random Genebody","Random GB-Enhancer"),
       col =  c("red","green","black","grey"),pch=16)
legend("topleft",legend =
       c(paste("n = ",sum(!is.na(cor.dist.rnd.gene$spearman))),
               paste("n = ",sum(!is.na(cor.dist.enh.rnd.gene$spearman)))),
       col = c("red","green"), pch = 16
       )
dev.off()
