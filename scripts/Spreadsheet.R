tx.red <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene"))
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes <- unlist(tx.red)
genes$gene_id <- names(genes)

prom.2kb <- promoters(genes,upstream=2000,downstream=0)
prom.tx.2kb <- promoters(tx,upstream=2000,downstream=0)
genebody.2kb <- promoters(genes,upstream=0,downstream=2000)

tss <- promoters(genes,upstream=1,downstream=0)

## SELECT MOST VARIABLE CpGs IN DIFFERENT GENOMIC COMPARTMENTS
selected.cg <- rownames(M.entity)

## promoter CpGs
ov <- findOverlaps(hm450[selected.cg],prom.2kb)
i.prom.cg <- unique(queryHits(ov))
selected.prom <- selected.cg[i.prom.cg]

## gene-body CpGs
ov <- findOverlaps(hm450[selected.cg],genebody.2kb)
i.gb.cg <- unique(queryHits(ov))
selected.gb <- selected.cg[i.gb.cg]

## distal CpGs
selected.dist <- selected.cg[-c(i.prom.cg,i.gb.cg)]



plot(density(cor.all$pearson,na.rm = T),xlim=c(-1,1),yaxt='n',ylab=" ")
lines(density(cor.prom$pearson,na.rm = T))
lines(density(cor.gb$pearson,na.rm = T),col=2)

lines(density(cor.gb1$pearson,na.rm = T),col=3)
lines(density(cor.promnew$pearson,na.rm = T),col=4)


#[which(cor.prom$spearman > 0.5),]


gr <- hm450[cor.prom$CG]
gr$mean <- STAT.entity[names(gr),]$mean
gr$cor <- cor.prom$spearman

gr <- gr[which(!values(gr)$cor %in% NA)]
values(gr)$cor
ccols  <- (as.numeric(as.factor(((gr$cor)))))

plotGrandLinear(gr, aes(y = mean, color =  ccols))
                

gr <- hm450[cor.gb$CG]
gr$mean <- STAT.entity[names(gr),]$mean
gr$cor <- cor.gb$spearman

gr <- gr[which(!values(gr)$cor %in% NA)]
values(gr)$cor
gr <- gr[seqnames(gr) %in% "chr1",]
ccols  <- (as.numeric(as.factor(((gr$cor)))))

plotGrandLinear(gr, aes(y = mean, color =  ccols)) + ggtitle("Correlation Genebody chr1")

library(ggbio) 
library(Homo.sapiens) 
class(Homo.sapiens) ##
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)
p.txdb <- autoplot(Homo.sapiens, which  = wh)
p.txdb
autoplot(Homo.sapiens, which  = wh, label.color = "black", color = "brown",
         fill = "brown")

tks <- tracks(gr, p.txdb,heights = c(2, 3)) + xlim(gr17) + theme_tracks_sunset()

tks

# gr <- hm450[cor.prom[which(cor.prom$spearman < -0.5),]$CG]
# gr$mean <- STAT.entity[names(gr),]$mean
# par(mfrow=c(2,1))
#plotGrandLinear(gr, aes(y = meanpos, meanneg), color = c("#7fc97f", "#fdc086"))


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggbio)
###
#Genebody
###
gr <- hm450[cor.gb$CG]
gr$mean <- STAT.entity[names(gr),]$mean
gr$correlation <- cor.gb$spearman

data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- transcriptsBy(txdb, by = "gene")
model17 <- subsetByOverlaps(model, genesymbol["HOXA5"]) 
exons <- gr
#exon17 <- subsetByOverlaps(exons, genesymbol["LY6G5C"]) 
exon17 <- exons[which(cor.gb$symb == "HOXA5")]
exon17 <- sort(exon17)

## reduce to make sure there is no overlap
## just for example
## suppose

p17 <- autoplot(txdb, genesymbol["HOXA5"])
plotRangesLinkedToData(exon17, stat.y = "correlation",stat.ylab = "Correlation",
                       annotation = list(p17), stat.coord.trans = coord_cartesian(),
                       theme.stat = theme_tracks_sunset(), theme.align = theme_clear(), 
                       linetype = 2)


genesymbol["DAXX"]


### ________________
###
#Promoter
###
gr <- hm450[cor.prom$CG]
gr$mean <- STAT.entity[names(gr),]$mean
gr$correlation <- cor.prom$spearman


data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- exonsBy(txdb, by = "tx")
model17 <- subsetByOverlaps(model, genesymbol["HOXA5"]) 
exons <- gr
#exon17 <- subsetByOverlaps(exons, genesymbol["TSPYL5"]) 
exon17 <- exons[which(cor.prom$symb == "HOXA5")]
## reduce to make sure there is no overlap
## just for example
## suppose
exon17 <- sort(exon17)
p17 <- autoplot(txdb, genesymbol["HOXA5"])
plotRangesLinkedToData(exon17, stat.y = "correlation",stat.ylab = "Correlation",
                       annotation = list(p17), stat.coord.trans = coord_cartesian(),
                       theme.stat = theme_tracks_sunset(), theme.align = theme_clear(), 
                       linetype = 2)

###
#Distal
###

gr <- (hm450[cor.dist$CG])

gr$mean <- STAT.entity[names(gr),]$mean
gr$correlation <- cor.dist$spearman

data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- transcriptsBy(txdb, by = "gene")
model17 <- subsetByOverlaps(model, genesymbol["HOXA5"]) 
exons <- gr
#exon17 <- subsetByOverlaps(exons, genesymbol["LY6G5C"]) 
exon17 <- exons[which(cor.dist$symb == "HOXA5")]
exon17 <- sort(exon17)
## reduce to make sure there is no overlap
## just for example
## suppose

p17 <- autoplot(txdb, genesymbol["HOXA5"],names.expr = "HOXA5",mode="reduce",geom = "alignment", ideogram = T )
plotRangesLinkedToData(exon17, stat.y = "correlation",stat.ylab = "Correlation",
                       annotation = list(p17), stat.coord.trans = coord_cartesian(),
                       theme.stat = theme_tracks_sunset(), theme.align = theme_clear(), 
                       linetype = 2)

####
#All Subsets, 3 Colors
####


cor.dist$type <- "dist"
cor.prom$type <- "prom"
cor.gb$type <- "gb"

nrow(cor.dist)
corr <- rbind(cor.dist,cor.gb,cor.prom)
nrow(corr)

corr <- corr[which(!is.na(corr$spearman)),]

gr <- (hm450[corr$CG])
  
gr$iqr <- STAT.entity[names(gr),]$iqr
#gr$type <- corr[corr$CG==names(gr),]$type
gr$dist <- 0
gr$prom <- 0
gr$gb <- 0
ndist <- nrow(corr[corr$type=="dist",])
ngb <- nrow(corr[corr$type=="gb",])
nprom <- nrow(corr[corr$type=="prom",])

gr[(ndist+1):(ndist+nprom)]$dist <- 0
gr[(ndist+nprom+1):(ndist+nprom+ngb)]$dist <- 0

gr[1:ndist]$prom <- 0
gr[(ndist+nprom+1):(ndist+nprom+ngb)]$prom <- 0

gr[1:ndist]$gb <- 0
gr[(ndist+1):(ndist+nprom)]$gb <- 0


gr[1:ndist]$dist <- corr[corr[corr$CG ==names(gr),]$type == "dist",]$spearman
gr[(ndist+1):(ndist+nprom)]$prom <- corr[corr[corr$CG ==names(gr),]$type == "prom",]$spearman
gr[(ndist+nprom+1):(ndist+nprom+ngb)]$gb <- corr[corr[corr$CG ==names(gr),]$type == "gb",]$spearman

#gr <- sort(gr)

data(genesymbol, package = "biovizBase")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
model <- transcriptsBy(txdb, by = "gene")
model17 <- subsetByOverlaps(model, genesymbol["HOXA5"]) 
exons <- gr
#exon17 <- subsetByOverlaps(exons, genesymbol["LY6G5C"]) 
exon17 <- exons[which(corr$symb == "HOXA5")]
exon17 <- sort(exon17)


# exon17[exon17$type=="prom"]$type <- "TRUE"
# exon17[exon17$type=="gb"]$type <- "FALSE"
# exon17[exon17$type=="dist"]$type <- NA

#exon17$type <- as.logical(as.factor(exon17$type))
## reduce to make sure there is no overlap
## just for example
## suppose

p17 <- autoplot(txdb, genesymbol["HOXA5"],names.expr = "HOXA5",mode="reduce",geom = "alignment", ideogram = T )
plotRangesLinkedToData(exon17, stat.y = c("dist","gb","prom","iqr"), size = 2,
                       annotation = list(p17),
                       #sig = "type", sig.col =  c("orange","lightblue"),
                       #stat.coord.trans = coord_cartesian(),
                       #theme.stat = theme_tracks_sunset(), theme.align = theme_clear(),
                       stat.align = "Test",
                       linetype = 2)



