#!/usr/bin/env Rscript
#===============================================================================
#' Author: Maria E. Calleja
#' Date: 2022/03/27
#' Trying cluster with:
#' salloc -c 6 --mem=40G -p medium --job-name=Cytopenias
#' srun -p medium --mem=40G --cpus-per-task=6 
#' module load R/4.0.0-foss-2018b
#' module load Seurat/4.0.0-foss-2018b-R-4.0.0
# 'R
#===============================================================================
## Packages and General Options
library(gwascat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(hues)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(org.Hs.eg.db)
options(stringsAsFactors = FALSE)
set.seed(777)
########################################################
# Global variables
########################################################
PARENT_DIR<-"/home/mecc/Downloads"
dir.create(file.path(PARENT_DIR,"Challenge"))
PROJECT_DIR<-file.path(PARENT_DIR,"Challenge")
dir.create(file.path(PROJECT_DIR,"figures"))
FIGURES_DIR<-file.path(PROJECT_DIR,"figures")
dir.create(file.path(PROJECT_DIR,"results"))
RESULTS_DIR<-file.path(PROJECT_DIR,"results")
#DATA_DIR<-file.path(PROJECT_DIR,"Data")
########################################################
# ## FUNCTIONS.
########################################################



#===============================================================================
## MAIN.
################################################################################
## 23andme RawData input
################################################################################
d <- read.table(file.path(PARENT_DIR,"example_genome.txt"),
                sep="\t", header=FALSE,
                colClasses=c("character", "character", "numeric", "character"),
                col.names=c("rsid", "chrom", "position", "genotype"))
tmp <- d$chrom
d$chrom = ordered(d$chrom, levels=c(seq(1, 22), "X", "Y", "MT"))
## It's never a bad idea to check your work
stopifnot(all(as.character(tmp) == as.character(d$chrom)))
ggplot(d) + geom_bar(aes(chrom))


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts(txdb)
tx.by.gene <- transcriptsBy(txdb, "gene")
tx.by.gene
columns(org.Hs.eg.db)
select(org.Hs.eg.db, keys="BRCA2", columns =c("ENTREZID", "SYMBOL", "GENENAME"), keytype="SYMBOL")
levels(d$chrom) <- paste("chr", c(1:22, "X", "Y", "M"), sep="")
my.snps <- with(d, GRanges(seqnames=chrom, 
                           IRanges(start=position, width=1), 
                           rsid=rsid, genotype=genotype)) 
apoe.i <- findOverlaps(tx.by.gene["675"], my.snps)
hits <- subjectHits(apoe.i) 
my.snps[hits]

if (length(grep("gwascat", search()))>0) detach("package:gwascat")
library(gwascat)
objects("package:gwascat")
data(ebicat37)
gwtrunc = ebicat37
topTraits(gwtrunc)
gwrngs.emd <- as.data.frame(elementMetadata(gwtrunc))
dm <- merge(d, gwrngs.emd, by.x="rsid", by.y="SNPS")

risk.alleles <- gsub("[^\\-]*-([ATCG?])", "\\1", dm$STRONGEST.SNP.RISK.ALLELE)

i.have.risk <- mapply(function(risk, mine) {
  risk %in% unlist(strsplit(mine, ""))}, risk.alleles, dm$genotype)

dm$i.have.risk <- i.have.risk
my.risk <- dm[dm$i.have.risk, ]
rel.cols <- c(colnames(d), "DISEASE.TRAIT", "RISK.ALLELE.FREQUENCY",
              "P.VALUE", "i.have.risk", "X95..CI..TEXT.")
head(my.risk[order(my.risk$RISK.ALLELE.FREQUENCY), rel.cols], 1)
dm[which(dm$rsid == "rs11880316"), "INITIAL.SAMPLE.DESCRIPTION"]
head(my.risk[grep("European", my.risk$Initial.Sample.Size), rel.cols], 30)

####VIZ
data(hg19IdeogramCyto, package = "biovizBase")
head(hg19IdeogramCyto)
autoplot(hg19IdeogramCyto, layout = "karyogram", cytobands = TRUE)
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
head(hg19)
hg19 <- autoplot(hg19, layout = "karyogram", cytobands = TRUE)


p <- plotOverview(hg19IdeogramCyto, cytoband=FALSE)
(elementMetadata(gwtrunc)$my.genotype <- 
    d$genotype[(match(elementMetadata(gwtrunc)$SNPS, d$rsid))])

elementMetadata(gwtrunc)$my.risk <- with(elementMetadata(gwtrunc), 
                                        mapply(function(risk, mine) {
                                          risk %in% unlist(strsplit(mine, ""))
                                        }, gsub("[^\\-]*-([ATCG?])", "\\1", "STRONGEST SNP-RISK ALLELE"), my.genotype))
elementMetadata(gwtrunc)
hg19 + geom_hotregion(gwtrunc, aes(color=my.risk))
