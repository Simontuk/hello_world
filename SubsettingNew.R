## Prequisities:
tx.red <- reduce(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene"))
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes <- unlist(tx.red)
genes$gene_id <- names(genes)

# Different Subsets:
# Promoter:
prom.2kb <- promoters(genes,upstream=2000,downstream=500)
#Promoter from transcripts
prom.tx.2kb <- promoters(tx,upstream=2000,downstream=0)
#Genebody as defined by start-region (why?)
#genebody.2kb <- promoters(genes,upstream=0,downstream=2000)

##Define TSS
tss <- promoters(genes,upstream=1,downstream=0)

## Promoter (2000 upstream, 500 downstream)
## SELECT MOST VARIABLE CpGs IN DIFFERENT GENOMIC COMPARTMENTS
selected.cg <- rownames(M.entity)

## promoter CpGs
ov <- findOverlaps(hm450[selected.cg],prom.2kb)
i.prom.cg <- unique(queryHits(ov))
selected.prom <- selected.cg[i.prom.cg]

## gene-body CpGs
ov <- findOverlaps(hm450[selected.cg],genes)
i.gb.cg <- unique(queryHits(ov))
selected.gb1 <- selected.cg[setdiff(i.gb.cg,i.prom.cg)]


# ov <- findOverlaps(hm450[selected.cg],tss)
# i.tss.cg <- unique(queryHits(ov))
# selected.tss <- selected.cg[i.tss.cg]

## distal CpGs
selected.dist <- selected.cg[-c(i.prom.cg,i.gb.cg)]


#cor.tss <- corMX(selected.tss,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL')
#cor.tss.rnd <- corMX(selected.tss,w=0,X.entity,M.entity,prom=prom.2kb,subset='ALL',randproc='gene')

