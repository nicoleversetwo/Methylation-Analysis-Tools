library(RColorBrewer)
library(gplots)
library(data.table)
library(StatePaintR)
library(funciVar)
library(GenomicRanges)
library(VariantAnnotation)
library(UpSetR)
norm_NRF2mut = fread("~/normal.NRF2mut.metilene_sorted.out", sep="\t")

#trying different approach here
#bg_NRF2WT_NRF2mut = fread("~/NRF2WT.NRF2mut.metilene_sorted.out", sep="\t")
bg_NRF2WT_NRF2mut = fread("~/NRF2WT.NRF2mut.metilene_sorted.version2.out", sep="\t")
fg_NRF2WT_NRF2mut <- fread("~/NRF2WT.NRF2mut.metilene.sorted.q.out_qval.0.05.out")

input_NRF2WT_NRF2mut = fread("~/metilene_NRF2mut_NRF2WT.input", sep="\t")

statehub.encode.aws <- "http://s3-us-west-2.amazonaws.com/statehub-trackhub/tracks/5813b67f46e0fb06b493ceb0/hg38/ENCODE/"
segmentation.files <- c(
  paste0(statehub.encode.aws,
         "neutrophil.8mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "mcf-7.16mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "k562.19mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "induced_pluripotent_stem_cell.7mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "hepatocyte.9mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "hela-s3.13mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "hct116.12mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "gm12878.12mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "fibroblast_of_lung.13mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "fibroblast_of_dermis.8mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "dohh2.8mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "cd14-positive_monocyte.9mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "cardiac_muscle_cell.9mark.segmentation.bed"),
  paste0(statehub.encode.aws,
         "bipolar_spindle_neuron.8mark.segmentation.bed"))
encode.segmentations <- GetSegmentations(files = segmentation.files)

esegs <- unlist(encode.segmentations)

esegs

mcols(esegs[esegs$state %in% c("EAR", "AR"), ])$state <- "Active Enhancer"
mcols(esegs[esegs$state %in% c("EPR", "EWR"), ])$state <- "Weak Enhancer"
mcols(esegs[esegs$state %in% c("PAR"), ])$state <- "Active Promoter"
mcols(esegs[esegs$state %in% c("PPR", "PPWR", "PWR"), ])$state <- "Weak Promoter"
mcols(esegs[esegs$state %in% c("CTCF"), ])$state <- "CTCF"
mcols(esegs[esegs$state %in% c("TRS"), ])$state <- "TRS"
mcols(esegs[esegs$state %in% c("HET"), ])$state <- "HET"
mcols(esegs[esegs$state %in% c("SCR"), ])$state <- "SCR"
mcols(esegs[esegs$state %in% c("PAR","PPR","PPWR", "PR", "PWR"), ])$state <- "Promoter"

gr_bg_NRF2WT_NRF2mut <- GRanges(seqnames = bg_NRF2WT_NRF2mut$V1,
                             ranges = IRanges(start = bg_NRF2WT_NRF2mut$V2,
                                              end = bg_NRF2WT_NRF2mut$V3),
                             strand = "*", Q = bg_NRF2WT_NRF2mut$V4, meanMeth = bg_NRF2WT_NRF2mut$V5, CpG = bg_NRF2WT_NRF2mut$V6 )

DMR_NRF2WT_NRF2mut[DMR_NRF2WT_NRF2mut$score>0,]

variants2 <- list(fg=gr_fg_NRF2WT_NRF2mut$score>0, bg=gr_bg_NRF2WT_NRF2mut)

NRFWT_NRFmut <- CalculateEnrichment(variants = list(fg = variants2$bg[variants2$bg$Q < 0.01 & variants2$bg$meanMeth >= 0.3, ], bg = variants2$bg), esegs, feature.type = "segmentations", CI = 0.8, return.overlaps = TRUE, strict.subset = TRUE)

PlotEnrichment(NRFWT_NRFmut$enrichment, block1 = "state", color.by = "sample", ncol = 12)
