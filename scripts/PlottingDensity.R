library(ggplot2)
require(grid)
######
#Plotting density / Histograms of COAD-data 
#####
cor.prom.spearman <- data.frame(cor.prom[,"spearman"])
cor.prom.spearman.gene <- data.frame(cor.prom.rnd.gene[,"spearman"])
cor.prom.spearman.sh <- data.frame(cor.prom.rnd.sh[,"spearman"])
colnames(cor.prom.spearman)[1] <- "spearman"
colnames(cor.prom.spearman.gene)[1] <- "spearman"
colnames(cor.prom.spearman.sh)[1] <- "spearman"


cor.prom.spearman$reg <- "Observed"
cor.prom.spearman.gene$reg <- "random gene"
cor.prom.spearman.sh$reg <- "shuffled gene"

head(cor.prom.spearman.gene)

cor.prom.plot <- rbind(cor.prom.spearman,cor.prom.spearman.gene,cor.prom.spearman.sh)

# head(cor.prom.plot)

# 
# hist(xx2,prob=T,xlim=c(1,-1),border=3)
# hist(xx,prob=T,add=T,border=1)
# hist(xx1,prob=T,add=T,border=2)
# 
# xxgg <- rbind(cor.prom.sort,cor.prom.rnd.gene.sort,cor.prom.rnd.sh.sort)


#GGPLOT
ggplot(cor.prom.plot,aes(spearman,fill=reg)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  ggtitle("Correlation Promoter Methylation") + 
  theme(legend.key.size = unit(2, "cm"),text =element_text(size=20))


####
#Genebody
####
cor.gb.spearman <- data.frame(cor.gb[,"spearman"])
cor.gb.spearman.gene <- data.frame(cor.gb.rnd.gene[,"spearman"])
cor.gb.spearman.sh <- data.frame(cor.gb.rnd.sh[,"spearman"])
colnames(cor.gb.spearman)[1] <- "spearman"
colnames(cor.gb.spearman.gene)[1] <- "spearman"
colnames(cor.gb.spearman.sh)[1] <- "spearman"


cor.gb.spearman$reg <- "Observed"
cor.gb.spearman.gene$reg <- "random gene"
cor.gb.spearman.sh$reg <- "shuffled gene"

#head(cor.gb.spearman.gene)

cor.gb.plot <- rbind(cor.gb.spearman,cor.gb.spearman.gene,cor.gb.spearman.sh)

ggplot(cor.gb.plot,aes(spearman,fill=reg)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  ggtitle("Correlation Genebody Methlyation") + 
  theme(legend.key.size = unit(2, "cm"),text =element_text(size=20))



####
#Distal
####
cor.dist.spearman <- data.frame(cor.dist[,"spearman"],stringsAsFactors = F)
cor.dist.spearman.gene <- data.frame(cor.dist.rnd.gene[,"spearman"],stringsAsFactors = F)
cor.dist.spearman.sh <- data.frame(cor.dist.rnd.sh[,"spearman"],stringsAsFactors = F)
colnames(cor.dist.spearman)[1] <- "spearman"
colnames(cor.dist.spearman.gene)[1] <- "spearman"
colnames(cor.dist.spearman.sh)[1] <- "spearman"


cor.dist.spearman$reg <- "Observed"
cor.dist.spearman.gene$reg <- "random gene"
cor.dist.spearman.sh$reg <- "shuffled gene"
cor.enh.spearman$reg <- "enhancer"

head(cor.dist.spearman.gene)
head(cor.dist.spearman.sh)

cor.dist.plot <- rbind(cor.dist.spearman,cor.dist.spearman.gene,cor.dist.spearman.sh,cor.enh.spearman)

tail(cor.dist.plot)

#cor.dist.plot$spearman <- as.integer(cor.dist.plot$spearman)
cor.dist.plot <- cor.dist.plot[!is.na(cor.dist.plot$spearman),]

ggplot(cor.dist.plot,aes(spearman,fill=reg)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  ggtitle("Correlation Distal Methylation") + 
  theme(legend.key.size = unit(2, "cm"),text =element_text(size=20))


####
#Enhancer
####
cor.enh.spearman <- data.frame(cor.enh[,"spearman"],stringsAsFactors = F)
cor.enh.spearman.gene <- data.frame(cor.enh.rnd.gene[,"spearman"],stringsAsFactors = F)
cor.enh.spearman.sh <- data.frame(cor.enh.rnd.sh[,"spearman"],stringsAsFactors = F)
colnames(cor.enh.spearman)[1] <- "spearman"
colnames(cor.enh.spearman.gene)[1] <- "spearman"
colnames(cor.enh.spearman.sh)[1] <- "spearman"


cor.enh.spearman$reg <- "Observed"
cor.enh.spearman.gene$reg <- "random gene"
cor.enh.spearman.sh$reg <- "shuffled gene"

head(cor.enh.spearman.gene)

cor.enh.plot <- rbind(cor.enh.spearman,cor.enh.spearman.gene,cor.enh.spearman.sh)

cor.enh.plot <- cor.enh.plot[!is.na(cor.enh.plot$spearman),]


ggplot(cor.enh.plot,aes(spearman,fill=reg)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  ggtitle("Correlation Enhancer Methylation") + 
  theme(legend.key.size = unit(2, "cm"),text =element_text(size=20))

cor.enh.spearman$reg <- "Enhancer"
cor.dist.spearman$reg <- "Distal"

cor.dist.enh <-  rbind(cor.enh.spearman,cor.dist.spearman)

ggplot(cor.dist.enh,aes(spearman,fill=reg)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + 
  ggtitle("Correlation Enhancer Methylation") + 
  theme(legend.key.size = unit(2, "cm"),text =element_text(size=20))


##GB and GB enhancers
par(mfrow=c(1,1))
hist(cor.gb.enh$spearman,prob=T,col=rgb(0,1,0,0.5))

hist(cor.gb$spearman,prob=T,col=rgb(1,0,0,0.5),add=T)

hist(cor.gb.rnd.gene$spearman,prob=T,col=rgb(0,0,1,0.5),add=T)

##Dist and Dist enhancers
par(mfrow=c(1,1))
hist(cor.dist.enh$spearman,prob=T,col=rgb(0,1,0,0.5))

hist(cor.dist$spearman,prob=T,col=rgb(1,0,0,0.5),add=T)

hist(cor.dist.rnd.gene$spearman,prob=T,col=rgb(0,0,1,0.5),add=T)
#cor.dist.spearman



cor.enh.spearman <- cor.enh.spearman[!is.na(cor.enh.spearman$spearman),]
lines(density(cor.enh.spearman$spearman))

head(cor.enh.spearman[1])
