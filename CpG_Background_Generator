library(BSgenome.Hsapiens.UCSC.hg38)
for(seqname in seqnames(BSgenome.Hsapiens.UCSC.hg38)) {
  query <- BSgenome.Hsapiens.UCSC.hg38[[seqname]]
  chr_matches <- matchPattern("CG", query)
  chr_matches <- reduce(chr_matches, min.gapwidth = 100) ## maximum gap between CpGs to be in the same window
  cg_counts <- vcountPattern("CG", chr_matches)
  chr_matches <- chr_matches[cg_counts >= 10] ## minimum number of CpGs to be considered a window
  if(length(chr_matches) < 1) next()
  chr_matches <- GRanges(seqname, ranges(chr_matches), nCpG = cg_counts[cg_counts >= 10], seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
  if(seqname == seqnames(BSgenome.Hsapiens.UCSC.hg38)[[1]]) {
    matches <- chr_matches
  } else {
    matches <- c(matches, chr_matches)
  }
}
