#' @title Histogram of annotations around boundaries
#'
#' @description This function take the output of `boundArea()` to produce a graph of the distribution of the genomic annotations around TAD boundaries (real distance).
#' `areahist()` measures the distances between the boundaries and the annotation features (start, end, closest or the center of annotations) and plot the distribution as an histogram.
#' `annot.boundary = "closest"`, use the closest border (start or stop) to measure the distances from the boundary independently to the annotation strand. Note that all annotations that overlap a boundary will have a distance of 0.
#'
#' @details Several features distances can be count. As an example (see examples for visualizations) of 2 genes which start at 20Kb of a TAD boundary (width of the 2 genes are 10Kb and are on the reverse an forward strand respectively).
#' -If `annot.boundary = "start"`, the distances from the boundary are 20Kb,
#' -If `annot.boundary = "end"`, the distances are 10 and 30Kb respectively,
#' -If `annot.boundary = "center"`, the distances are 15Kb and 25Kb respectively,
#' -If `annot.boundary = "closest"`, the distances are 10Kb and 20Kb respectively.
#'
#' @inheritParams areaCov
#' @param annot.boundary Type of feature to analyzed. `"Start"`, `"end"` or `"center"` of each `annot.gr`. It is possible to use `"closest"` to take the closest border (start or end) independently to the strands.
#' @param annot.strand Logical. Default is `FALSE` to plot the distribution as histogram. If `TRUE`, distributions are separated according to their strands and are display with lines.
#' @param bin.width Size of the bin in base pair to count the number of annotations features. It should match the bin.width of the matrix used to call the domains. Default is `NULL` to to use a size in 10 time smaller than `window.size` parameter of `TADarea()`.
#'
#' @return `ggplot`
#'
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @importFrom IRanges distance IRanges ranges
#' @import ggplot2
#' @importFrom scales unit_format
#'
#' @export
#'
#' @examples
#' # create GRange with 2 annotations (genes with same size and start but different strand)
#' annot.gr <- tad.gr <- dataframes2grange(
#'   data.frame(chr = 1, start = c(20e3, 30e3), end = c(30e3, 40e3), strand = c("-", "+")),
#'   data.frame(chr = "1", size = 60e3),
#'   strand.col = 4
#' )
#'
#' # create GRange with 1 TAD
#' tad.gr <- dataframes2grange(
#'   data.frame(chr = 1, start = c(10e3), end = c(50e3)),
#'   data.frame(chr = "1", size = 60e3)
#' )
#'
#' # Vizualisation
#' TADplot(tad.gr = tad.gr, start = 5e3, stop = 55e3, chr = 1, annot.gr = annot.gr)
#'
#' # Distribution analysis surrounding the start of the TAD
#' data.gr <- boundArea(
#'   domain.gr = tad.gr,
#'   annot.gr = annot.gr,
#'   window.size = 50e3,
#'   domain.boundary = "start"
#' )
#'
#' # Visualized differences according to annotation features:
#' for (features in c("start", "end", "center", "closest")) {
#'   output <- areaHist(
#'     data.gr = data.gr,
#'     annot.strand = FALSE,
#'     bin.width = 1e3,
#'     annot.boundary = features
#'   )
#'   print(output)
#' }
areaHist <- function(data.gr, annot.boundary = "center", annot.strand = TRUE, bin.width = NULL) {

  window.size <- unique(data.gr$window.size) # window.size used in TADarea function
  GenomeInfoDb::seqlengths(data.gr) <- NA # remove seqlengths to prevent out-of-bound ranges

  if (is.null(bin.width)) {
    bin <- window.size / 10
  } else {
    bin <- bin.width # binwidth for ggplot
  }

  if (annot.boundary == "closest") {
    data2.gr <- data.gr
  }

  if (annot.boundary == "start" | annot.boundary == "end" | annot.boundary == "center") {
    data2.gr <- GenomicRanges::resize(data.gr, 1, fix = annot.boundary) # resize annotations according to start, stop or center

    data2.gr <- GenomicRanges::restrict(data2.gr, start = 1, end = 2 * window.size - 1) # remove ranges outside of the window
  }

  ###################################################################
  # mesure distances from annotation features from TAD borders
  ###################################################################
  data2.gr$distance <- ifelse(
    GenomicRanges::start(data2.gr) < window.size, # if start before TAD border
    -IRanges::distance(IRanges::IRanges(start = window.size), IRanges::ranges(data2.gr)), # negatives distances
    IRanges::distance(IRanges::IRanges(start = window.size), IRanges::ranges(data2.gr)) # else positive distances
  )

  ###################################################################
  # histo
  ###################################################################
  names(data2.gr) <- NULL # some names are duplicated, so remove the names to create a data.frame

  if (isFALSE(annot.strand)) {
    p <- ggplot2::ggplot(as.data.frame(data2.gr), ggplot2::aes(x = distance)) +
      ggplot2::geom_histogram(binwidth = bin, boundary = bin / 2, size = 0.75, color = "black", fill = "white")
  } # histogram

  if (isTRUE(annot.strand)) {
    p <- ggplot2::ggplot(as.data.frame(data2.gr), ggplot2::aes(x = distance, color = strand)) +
      ggplot2::geom_freqpoly(binwidth = bin, boundary = bin / 2, size = 0.75, alpha = 0.75) +
      ggplot2::scale_color_manual(values = c("#33A02C", "#1F78B4")) +
      ggplot2::geom_point(stat = "bin", aes(y = after_stat(count)), binwidth = bin, boundary = bin / 2, size = 1, alpha = 0.75)
  }

  return(
    p + ggplot2::scale_x_continuous(
      labels = scales::unit_format(unit = "kb", scale = 1e-3),
      breaks = c(
        seq(-window.size + bin / 2, -bin / 2, 1 * 10^(nchar(window.size - 1) - 1) * 2), 0,
        rev(seq(-window.size + bin / 2, -bin / 2, 1 * 10^(nchar(window.size - 1) - 1) * 2)) * -1
      ), # adjust the x-axis according to the size of the bins
      limits = c(-window.size + bin / 2, window.size - bin / 2)
    ) +
      ggplot2::geom_vline(aes(xintercept = 0), color = "red", size = 0.75, alpha = 0.75) +
      ggtitle(annot.boundary)
  )
}
