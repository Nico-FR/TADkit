#' @title Compartments calling (A or B) from PC1 scores
#'
#' @description From the scores of each bin (i.e first principal component scores), `PC1calling` output a `GRanges` with compartment A for score >= 0 and compartment B for score < 0
#'
#' @details From bedgraph as input (`GRanges` object or path) output a `GRanges` with the domains A or B.
#' `NA` score are considered as part of a compartment A if they are between two compartments A (the other way around for B compartments).
#' Otherwise (if `NA` are between compartment A and B) they are not called, thus leaving a gap between two compartments.
#'
#'
#' @param bedgraph `GRanges` or `data.frame` object or bedgraph path (data frame without header and 4 columns tab separated).
#'
#' @return `GRanges`
#'
#' @import GenomeInfoDb
#' @importFrom BiocGenerics sort
#' @importFrom dplyr mutate case_when
#' @importFrom S4Vectors split
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' PC1calling(PC1_250kb.gr)
#'

PC1calling <- function(bedgraph) {
  if (is.character(bedgraph)) {
    bedgraph = read.table(bedgraph, header = FALSE, sep = "\t") %>% dplyr::mutate(comp = dplyr::case_when(V4 < 0 ~ 'B', V4 >= 0 ~ 'A', is.na(V4) ~ "AB"))

    grange = GenomicRanges::makeGRangesFromDataFrame(bedgraph, start.field="V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
  }

  if (is.data.frame(bedgraph)) {
    grange = bedgraph %>% dplyr::mutate(comp = dplyr::case_when(bedgraph[,4] < 0 ~ 'B', bedgraph[,4] >= 0 ~ 'A', is.na(bedgraph[,4]) ~ "AB")) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  }

  seqinfo <- NULL
  if (inherits(bedgraph, "GRanges")) {
    seqinfo = as.data.frame(GenomeInfoDb::seqinfo(bedgraph))
    grange = bedgraph %>% as.data.frame %>% dplyr::mutate(comp = dplyr::case_when(GenomicRanges::mcols(bedgraph)[,1] < 0 ~ 'B', GenomicRanges::mcols(bedgraph)[,1] >= 0 ~ 'A', is.na(GenomicRanges::mcols(bedgraph)[,1]) ~ "AB")) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  }



  gap = GenomicRanges::gaps(grange, start = 2)
  if (length(gap) > 0) {
    gap$comp = "AB"
    grange = c(grange, gap + 1) %>% sort()
  }


  data1 = S4Vectors::split(grange, grange$comp)

  if (length(data1$AB) > 0){
    B = GenomicRanges::reduce(c(data1$B, GenomicRanges::reduce(data1$AB)[GenomicRanges::countOverlaps(GenomicRanges::reduce(data1$AB), data1$B, maxgap = 1) == 2]))

    A = GenomicRanges::reduce(c(data1$A, GenomicRanges::reduce(data1$AB)[GenomicRanges::countOverlaps(GenomicRanges::reduce(data1$AB), data1$A, maxgap = 1) == 2]))
  }

  if (length(data1$AB) == 0){
    B = GenomicRanges::reduce(data1$B)

    A = GenomicRanges::reduce(data1$A)
  }


  A$comp = "A"
  B$comp = "B"
  output.gr = BiocGenerics::sort(c(A,B))
  names(output.gr) = output.gr$comp

  #add seqinfo if available
  if (!is.null(seqinfo)) {
    for (i in 1:length(GenomeInfoDb::seqlengths(output.gr))) {
      GenomeInfoDb::seqlengths(output.gr)[i] <- as.numeric(seqinfo[, 1][rownames(seqinfo) == names(GenomeInfoDb::seqlengths(output.gr))[i]])
    }}

  return(output.gr)
}
