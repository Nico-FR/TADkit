#' Convert boundary to domain
#'
#' @description From datas with boundaries stored by line (chromosome, start and end), `boundary2domain()` return `GRanges` with domains.
#'
#' @details Start and end of domains are the middle of boundaries.
#'
#' @param boundaries `dataframe` (chromosome, start and end) or `GRanges` with the boundaries.
#'
#' @return `GRanges` object with domains
#'
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' #create GRanges with chr sizes:
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#'
#' domains.gr = boundary2domain(boundaries.gr)
#' domains.gr
#'
boundary2domain <- function(boundaries) {

  seqinfo <- NULL
  if (inherits(boundaries, "GRanges")) {
    seqinfo = as.data.frame(GenomeInfoDb::seqinfo(boundaries))
    boundaries = as.data.frame(boundaries)[,1:3]}

  boundaries$mean = apply(boundaries[,2:3], 1, mean)
  list = split(boundaries, boundaries[,1])

  df.lst = lapply(list, function(x) {
    data.frame(
      seqnames = x[1:(nrow(x) - 1), 1],
      start = x$mean[1:(nrow(x) - 1)],
      end = x$mean[2:nrow(x)])
    })
  df = do.call("rbind", df.lst)
  rownames(df) = paste0(df$seqnames, "_", df$start)
  gr = GenomicRanges::makeGRangesFromDataFrame(df)

  #add seqinfo if available
  if (!is.null(seqinfo)) {
    for (i in 1:length(GenomeInfoDb::seqlengths(gr))) {
      GenomeInfoDb::seqlengths(gr)[i] <- as.numeric(seqinfo[, 1][rownames(seqinfo) == names(GenomeInfoDb::seqlengths(gr))[i]])
    }}
  return(gr)
}
