# cormax = 1
# cormin = 0.5
CG.GENE <- cor.dist#[which(cor.dist$spearman < cormax & cor.dist$spearman > cormin),]
CG.gr <- hm450[CG.GENE$CG]
length(CG.gr)
class(CG.gr)
load("/icgc/dkfzlsdf/analysis/B080/steiger/data/roadmap/roadmap.seg.enh.RData")
class(roadmap.seg.enh)
library(GenomicRanges)

roadmap.seg.enh <- GRangesList(roadmap.seg.enh)
CGenh.gr <- subsetByOverlaps(CG.gr,roadmap.seg.enh)

cor.enh <- cor.dist[cor.dist$CG %in% names(CGenh.gr),]
cor.enh.rnd.gene <- cor.dist.rnd.gene[cor.dist.rnd.gene$CG %in% names(CGenh.gr),]
cor.enh.rnd.sh <- cor.dist.rnd.sh[cor.dist.rnd.sh$CG %in% names(CGenh.gr),]

dim(cor.enh)
