#------------------------------------------------------------
# Auxiliary function
#------------------------------------------------------------
# o input: file with region
# o output: list with region dataframe, region genomicRanges
# o description: 
#   1. read a DMR region
#   2. map to 10 upstream and 10 downstream genes
#   3. annotate genes to promoter (<2kb from region) or distal (>=2kb from region)
#   4. Add nearest gene to 
read_region <- function(file){
  region <- readr::read_tsv(file,col_names = FALSE,col_types = readr::cols())
  colnames(region) <- c("chr","start","end","score")
  
  region$ID <- paste0(region$chr,"_",
                      region$start,"_",
                      region$end)
  
  region.gr <- region %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepStandardChromosomes(pruning.mode = 'coarse')
  
  names(region.gr) <- region.gr$ID
  
  gene.info <- TCGAbiolinks:::get.GRCh.bioMart(genome = "hg38",as.granges = TRUE)
  # Get ensemble gene information (at transcript level) - used to calculate distance
  tssAnnot <- ELMER::getTSS(genome = "hg38")
  # Use ELMER function
  suppressWarnings({
    near.genes.table <- ELMER:::getRegionNearGenes(TRange = region.gr,
                                                   geneAnnot = gene.info,
                                                   numFlankingGenes = 20,
                                                   tssAnnot = tssAnnot)
    
  })
  
  # Annotate region as distal (no tss close to 2kb) or promoter (tss close to 2kb)
  promoter <- promoters(tssAnnot,upstream = 2000,downstream = 2000)
  region.gr$Away_2KB_TSS <- "Distal (> 2KB away TSS)"
  region.gr$Away_2KB_TSS[region.gr$ID %in% names(subsetByOverlaps(region.gr,promoter))] <- "Promoter (< 2KB away TSS)"
  
  # Add distance from region to the nearest tss
  tss.start <- promoters(tssAnnot,upstream = 1,downstream = 1)
  region.gr$distance_nearest_tss <- distanceToNearest(region.gr,tss.start) %>% values %>%  as.data.frame() %>% pull()
  
  # get neastest gene information
  values(region.gr) <- cbind(values(region.gr),tssAnnot[nearest(region.gr,tssAnnot)] %>% values)
  
  # output:
  # 1) the region read with neastest gene information as GRanges
  # 2) a table mapping each region to 10 upstream and 10 downstream genes
  return(list("region" = region.gr,
              "near.genes.table" = near.genes.table))
}
