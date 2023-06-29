##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
##   Import Files
##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
​
SA09367.meth <- read.delim("~/SA09367.sorted.bed.hg38_112PMD", 
                                 header=FALSE, 
                                 stringsAsFactors=FALSE)
​
SA09365.meth <- read.delim("~/SA09365.sorted.bed.hg38_112PMD", 
                           header=FALSE, 
                           stringsAsFactors=FALSE)
​
​
SA09367.read <- read.delim("~/SA09367.sorted.coverage.bed.hg38_112PMD", 
                                          header=FALSE,
                                          stringsAsFactors = FALSE)
​
SA09365.read <- read.delim("~/SA09365.sorted.coverage.bed.hg38_112PMD", 
                           header=FALSE,
                           stringsAsFactors = FALSE)
​
##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
##   Manipulate Files for MethylSeekR
##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
​
library(plyr)
library(dplyr)
library(tidyr)
​
SA09365.df <- mutate(SA09365.read, pos = SA09365.read$V3 - 1)
SA09365.df <- SA09365.df[,c(1,5,4)]
SA09365.df <- mutate(SA09365.df, count = SA09365.read$V4 * SA09365.meth$V4 )
SA09365.df <- as.data.frame(SA09365.df)
colnames(SA09365.df) <- c("chr","pos","T","M")
#SA09365.df <- complete.cases(SA09365.df)
SA09365.df$T <- as.numeric(SA09365.df$T)
​
​
SA09365.df$pos <- as.integer(SA09365.df$pos)
SA09365.df$V1 <- as.character(SA09365.df$V1)
​
SA09367.df <- mutate(SA09367.read, pos = SA09367.read$V3 - 1)
SA09367.df <- SA09367.df[,c(1,5,4)]
SA09367.df <- mutate(SA09367.df, count = SA09367.read$V4 * SA09367.meth$V4 )
colnames(SA09367.df) <- c("chr","pos","T","M")
#SA09367.df <- complete.cases(SA09367.df)
​
##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
##   MethylSeekR
##~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
​
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg38")
library(MethylSeekR)
set.seed(123)
​
sLengths=seqlengths(Hsapiens)
head(sLengths)
library(GenomicRanges)
meth.9365.gr <- makeGRangesFromDataFrame(df2,
                                 start.field = "pos",
                                 end.field = "pos",
                                 seqnames.field = "chr",
                                 keep.extra.columns = T)
​
PMDsegments.gr <- segmentPMDs(m=meth.9365.gr, 
                              chr.sel="chr16",
                              seqLengths=sLengths, 
                              num.cores=1,
                              pdfFilename = "~/PMD_chr16-2.pdf")
​
plotAlphaDistributionOneChr(m=meth.9365.gr, 
                            chr.sel="chr16",
                            num.cores=1, 
                            pdfFilename = "~/AlphaDistribution_chr16.pdf")
​
plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.9365.gr,
                                               PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), 
                            chr.sel="chr16",
                            num.cores=1,
                            pdfFilename = "~/chr16_notPMD_seq_7.25.19.pdf")
​
​
​
​
plotAlphaDistributionOneChr(m=subsetByOverlaps(meth.9365.gr,
                                               PMDsegments.gr[values(PMDsegments.gr)$type=="notPMD"]), 
                            chr.sel="chr16",
                            num.cores=1)
​
plotPMDSegmentation(m=meth.9365.gr, segs=PMDsegments.gr,
                    pdfFilename = "~/chr16_PMD_seq_7.25.19.pdf")
​
savePMDSegments(PMDs=PMDsegments.gr,
                  GRangesFilename="PMDs.gr.rds", 
                  TableFilename="PMDs.tab")
