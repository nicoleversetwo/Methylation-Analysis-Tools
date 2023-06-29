library(GenomicRanges)
library(tidyr)
library(dplyr)
# get beta values
beta <- readr::read_tsv("~/Downloads/All_M_8.21.18.txt")
beta <- All_M_8.21.18

cpgs <- makeGRangesFromDataFrame(beta,
                                 start.field = "ESCC_9_M",
                                 end.field = "ESCC_9_M",
                                 seqnames.field = "ESCC_10_M",
                                 keep.extra.columns = T)

# get regions
dmrseq <- readr::read_tsv("~/Downloads/dmrseq_Sig_Hypo_ESCC__ESCCvsEAC_hg38_8.7.18.txt",
                          col_names = TRUE)
dmrseq <- merged
colnames(dmrseq) <- c("ID",colnames(dmrseq)[-14])
regions <- makeGRangesFromDataFrame(dmrseq)

# for each region identify probes inside it and caculate the avarage
ret <- adply(1:length(regions),
             .margins = 1,
             .fun = function(x){
               # get cpgs inside region and calculte mean by samples
               overlapping.cpgs <- as.matrix(values(subsetByOverlaps(cpgs,regions[x])))
               colMeans(overlapping.cpgs,na.rm = T) 
             },.progress = "time")
ret$regions <- as.data.frame(regions) %>% unite("ID","seqnames","start","end") %>% pull(ID)
