#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh.bioMart <- function(genome="hg19") {
  tries <- 0L
  msg <- character()
  while (tries < 3L) {
    gene.location <- tryCatch({
      if (genome == "hg19"){
        # for hg19
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           host = "feb2014.archive.ensembl.org",
                           path = "/biomart/martservice" ,
                           dataset = "hsapiens_gene_ensembl")
        attributes <- c("chromosome_name",
                        "start_position",
                        "end_position", "strand",
                        "ensembl_gene_id", "entrezgene",
                        "external_gene_id")
      } else {
        # for hg38
        ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        attributes <- c("chromosome_name",
                        "start_position",
                        "end_position", "strand",
                        "ensembl_gene_id",
                        "external_gene_name")
      }
      description <- listDatasets(ensembl)[listDatasets(ensembl)$dataset=="hsapiens_gene_ensembl",]$description
      message(paste0("Downloading genome information (try:", tries,") Using: ", description))
      
      filename <-  paste0(gsub("[[:punct:]]| ", "_",description),".rda")
      if(!file.exists(filename)) {
        chrom <- c(1:22, "X", "Y")
        gene.location <- getBM(attributes = attributes,
                               filters = c("chromosome_name"),
                               values = list(chrom), mart = ensembl)
        save(gene.location, file = filename)
      } else {
        message("Loading from disk")
        gene.location <- get(load(filename))
      }
      gene.location
    }, error = function(e) {
      msg <<- conditionMessage(e)
      tries <<- tries + 1L
    })
    if(!is.null(gene.location)) break
  }
  if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
  
  return(gene.location)
}
