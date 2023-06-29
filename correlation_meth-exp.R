# Update ELMER package
# devtools::install_github("tiagochst/ELMER")
library(ELMER) 
library(GenomicRanges)
​
#------------------------------------------------------------
# Get nearest genes HYPO in EAC
#------------------------------------------------------------
​
dmr <- (Hypofiltered)
​
dmr$ID <- paste0(dmr$seqnames,":",
                 dmr$start,"-",
                 dmr$end)
​
​
dmr.gr <- makeGRangesFromDataFrame(dmr,keep.extra.columns = T)
​
# Get ensemble gene information (at gene level)
gene.info <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38",
                                             as.granges = TRUE)
​
names(dmr.gr) <- dmr.gr$ID
​
# How it works: for each region in TRange the n nearest regions (from geneAnnoted) will be mapped
# output: list
# input: two granges.
dmr.near.genes <- ELMER::GetNearGenes(TRange = dmr.gr, 
                                      geneAnnot = gene.info,
                                      numFlankingGenes = 20,
                                      cores = 1)
​
# Make list to data frame
dmr.near.genes.table <- data.table::rbindlist(dmr.near.genes)
​
write.table(dmr.near.genes.table, file = "ALL_HYPO_nearGenesTable_7.2.18.txt",append = FALSE, quote = FALSE, sep = "\t",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
​
#------------------------------------------------------------
# Get nearest genes HYPO in ESCC
#------------------------------------------------------------
​
dmr2 <- (Hyperfiltered)
​
dmr2$ID <- paste0(dmr2$seqnames,":",
                 dmr2$start,"-",
                 dmr2$end)
​
​
dmr2.gr <- makeGRangesFromDataFrame(dmr2,keep.extra.columns = T)
​
# Get ensemble gene information (at gene level)
gene.info <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38",
                                             as.granges = TRUE)
​
names(dmr2.gr) <- dmr2.gr$ID
​
# How it works: for each region in TRange the n nearest regions (from geneAnnoted) will be mapped
# output: list
# input: two granges.
dmr2.near.genes <- ELMER::GetNearGenes(TRange = dmr2.gr, 
                                      geneAnnot = gene.info,
                                      numFlankingGenes = 20,
                                      cores = 1)
​
# Make list to data frame
dmr2.near.genes.table <- data.table::rbindlist(dmr2.near.genes)
​
write.table(dmr2.near.genes.table, file = "ALL_HYPER_nearGenesTable_7.2.18.txt",append = FALSE, quote = FALSE, sep = "\t",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
​
#------------------------------------------------------------
# Correlated expression with t.test for each gene, Hypo in EAC
#------------------------------------------------------------
# donwload gene expression
library(dplyr)
library(plyr)
library(SummarizedExperiment)
ELMER::getTCGA("ESCA",RNA = TRUE,Meth = FALSE)
load("data/ESCA/ESCA_RNA_hg38_no_filter.rda")
​
esad.id <- readr::read_csv("~/esad_id.csv") %>% pull(2)
escc.id <- readr::read_csv("~/escc_id.csv") %>% pull(2)
data$group <- NA
data$group[colnames(data) %in% esad.id] <- "ESAD"
data$group[colnames(data) %in% escc.id] <- "ESCC"
​
exp <- assay(data)
methy <- which(data$group == "ESCC")
unmethy <- which(data$group == "ESAD")
​
# for each gene we will test if it is up regulated in the hypo methylated group
Probe.gene <- adply(.data = unique(dmr.near.genes.table$GeneID), .margins = 1,
                    .fun = function(gene) {
                      if(!gene %in% rownames(exp)) return(NA)
                      # is the exp in the unmethylated group greater than methylated ?
                      t.test(x = exp[gene,unmethy],y = exp[gene,methy],alternative = "greater")$p.value
                    },
                    .progress = "text", .parallel = FALSE, .id = NULL
)
colnames(Probe.gene) <- "exp_pval"
Probe.gene$exp_FDR <- p.adjust(Probe.gene$exp_pval,method = "BH")
Probe.gene$GeneID <- unique(dmr.near.genes.table$GeneID)
Probe.gene$status <- "Insignificant"
Probe.gene$status[Probe.gene$exp_FDR < 0.01] <- "Negative correlated"
​
correlation <- merge(dmr.near.genes.table,Probe.gene,by = "GeneID")
​
colnames(correlation)[2] <- "ID"
​
# we wanna have a single object with expression t.test and region methylation analysis
final <- merge(dmr,correlation,by = "ID")
​
write.table(final, file = "HypoInEAC_7.3.18.txt",append = FALSE, quote = FALSE, sep = "\t",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
​
#------------------------------------------------------------
# Correlated expression with t.test for each gene, Hypo in ESCC
#------------------------------------------------------------
# donwload gene expression
library(dplyr)
library(plyr)
library(SummarizedExperiment)
ELMER::getTCGA("ESCA",RNA = TRUE,Meth = FALSE)
load("data/ESCA/ESCA_RNA_hg38_no_filter.rda")
​
esad.id <- readr::read_csv("~/esad_id.csv") %>% pull(2)
escc.id <- readr::read_csv("~/escc_id.csv") %>% pull(2)
data$group <- NA
data$group[colnames(data) %in% esad.id] <- "ESAD"
data$group[colnames(data) %in% escc.id] <- "ESCC"
​
exp <- assay(data)
methy <- which(data$group == "ESAD")
unmethy <- which(data$group == "ESCC")
​
# for each gene we will test if it is up regulated in the hypo methylated group
Probe.gene <- adply(.data = unique(dmr2.near.genes.table$GeneID), .margins = 1,
                    .fun = function(gene) {
                      if(!gene %in% rownames(exp)) return(NA)
                      # is the exp in the unmethylated group greater than methylated ?
                      t.test(x = exp[gene,unmethy],y = exp[gene,methy],alternative = "greater")$p.value
                    },
                    .progress = "text", .parallel = FALSE, .id = NULL
)
colnames(Probe.gene) <- "exp_pval"
Probe.gene$exp_FDR <- p.adjust(Probe.gene$exp_pval,method = "BH")
Probe.gene$GeneID <- unique(dmr2.near.genes.table$GeneID)
Probe.gene$status <- "Insignificant"
Probe.gene$status[Probe.gene$exp_FDR < 0.01] <- "Negative correlated"
​
correlation <- merge(dmr2.near.genes.table,Probe.gene,by = "GeneID")
​
colnames(correlation)[2] <- "ID"
​
# we wanna have a single object with expression t.test and region methylation analysis
final2 <- merge(dmr2,correlation,by = "ID")
​
write.table(final2, file = "HypoInESCC_7.3.18.txt.txt",append = FALSE, quote = FALSE, sep = "\t",
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
