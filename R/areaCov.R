#' @title Density plot of annotations around boundaries
#'
#' @description This function take the output of `boundArea()` to produce a graph of the density of the genomic annotations around TAD boundaries (real distance).
#'
#' @details If the density is significantly different between strands, it can be useful to use relative density (i.e normalized using the zscore with norm = TRUE parameter).
#'
#' @param data.gr Output of `boundArea()`.
#' @param annot.strand Logical. If `TRUE` (default): separate coverage according to their strands.
#' @param bin.width Size of the sliding window in base pair to calculate the mean coverage. Default is `NULL` to use a size 10 times smaller than `window.size` parameter of `boundArea()`.
#' @param norm Logical. Normalized density using relative content (zscore of density). Default is `FALSE`.
#'
#' @return `ggplot`
#'
#' @importFrom S4Vectors runmean
#' @importFrom IRanges coverage ranges
#' @import GenomicRanges
#' @import ggplot2
#' @importFrom scales unit_format
#'
#' @export
#'
#' @examples
#' # create GRange with 1 TAD
#' tad.gr <- dataframes2grange(
#'   data.frame(chr = 1, start = c(10e3), end = c(50e3)),
#'   data.frame(chr = "1", size = 60e3)
#' )
#'
#' # create annotations of 12 genes every 2kb with width = 10kb on both strand:
#' annot.gr <- dataframes2grange(
#'   data.frame(
#'     chr = rep("1", 12),
#'     start = seq(20e3, 42e3, 2e3),
#'     end = seq(20e3, 42e3, 2e3) + 10e3,
#'     strand = c(rep("+", 6), rep("-", 6))
#'   ),
#'   data.frame(chr = "1", size = 60e3),
#'   strand.col = 4
#' )
#' # Vizualisation
#' TADplot(tad.gr = tad.gr, start = 5e3, stop = 55e3, chr = 1, annot.gr = annot.gr)
#'
#' # Distribution analysis surrounding the start of the TAD
#' data.gr <- boundArea(
#'   domain.gr = tad.gr,
#'   annot.gr = annot.gr,
#'   window.size = 60e3,
#'   domain.boundary = "start"
#' )
#'
#' # plot coverage
#' areaCov(data.gr = data.gr, annot.strand = FALSE)
#' areaCov(data.gr = data.gr, annot.strand = TRUE)
#'
areaCov <- function(data.gr, annot.strand = TRUE, bin.width = NULL, norm = FALSE) {
  window.size <- unique(data.gr$window.size) # window.size used in TADarea function

  if (is.null(bin.width)) {
    bin <- unique(data.gr$window.size) / 10
  } else {
    bin <- bin.width # binwidth for ggplot
  }

  # cut GRanges spanning the window
  data.gr <- GenomicRanges::restrict(data.gr, start = 1, end = 2 * window.size - 1)

  #####################################################################################
  # coverage on sliding windows (strand + & - separated)
  #####################################################################################
  if (isTRUE(annot.strand)) {
    # coverage for + strand
    cov.rle <- IRanges::coverage(IRanges::ranges(data.gr[BiocGenerics::strand(data.gr) == "+"]),
      width = 2 * window.size
    )
    mean_cov.rle <- S4Vectors::runmean(cov.rle, bin) # coverage for sliding window of length = bin
    mean_cov.table.f <- as.data.frame(mean_cov.rle) # translate to data.frame
    mean_cov.table.f$index <- c(1:(length(mean_cov.table.f$value)) - window.size + bin / 2)
    mean_cov.table.f$strand <- rep("+", length(mean_cov.table.f$value))
    mean_cov.table.f$zscore <- scale(mean_cov.table.f$value) # normalization using zscore
    mean_cov.table.f$density <- mean_cov.table.f$value / data.gr$nb_boundary[1] # transforms coverage into density

    # coverage for - strand
    cov.rle <- IRanges::coverage(IRanges::ranges(data.gr[BiocGenerics::strand(data.gr) == "-"]),
      width = 2 * window.size
    )
    mean_cov.rle <- S4Vectors::runmean(cov.rle, bin) # coverage for sliding window of 10kb
    mean_cov.table.r <- as.data.frame(mean_cov.rle) # translate to data.frame
    mean_cov.table.r$index <- c(1:(length(mean_cov.table.r$value)) - window.size + bin / 2)
    mean_cov.table.r$strand <- rep("-", length(mean_cov.table.r$value))
    mean_cov.table.r$zscore <- scale(mean_cov.table.r$value) # normalization using zscore
    mean_cov.table.r$density <- mean_cov.table.r$value / data.gr$nb_boundary[1] # transforms coverage into density

    # merge strands
    mean_cov.table <- rbind(mean_cov.table.r, mean_cov.table.f)
  }
  #####################################################################################
  # coverage on sliding windows (strand + & - merged)
  #####################################################################################
  if (isFALSE(annot.strand)) {
    cov.rle <- IRanges::coverage(IRanges::ranges(data.gr), width = 2 * window.size)
    mean_cov.rle <- S4Vectors::runmean(cov.rle, bin) # coverage for sliding window of 10kb
    mean_cov.table <- as.data.frame(mean_cov.rle) # translate to data.frame
    mean_cov.table$index <- c(1:(length(mean_cov.table$value)) - window.size + bin / 2)
    mean_cov.table$strand <- rep("*", length(mean_cov.table$value))
    mean_cov.table$zscore <- scale(mean_cov.table$value)
    mean_cov.table$density <- mean_cov.table$value / data.gr$nb_boundary[1] # transforms coverage into density
  }

  data <- mean_cov.table # input for ggplot
  names(data) <- c("coverage", "distance", "strand", "zscore", "density")

  # y values if norm = TRUE or FALSE
  if (isFALSE(norm)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(y = density, x = distance, color = strand))
  }
  if (isTRUE(norm)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(y = zscore, x = distance, color = strand))
  }

  # ggplot to return
  return(p + ggplot2::geom_vline(ggplot2::aes(xintercept = 0), color = "red", size = 0.75, alpha = 0.75) +
    ggplot2::geom_line(alpha = 0.5, size = 1) +
    ggplot2::xlab(data.gr$boundary_type[1]) +
    ggplot2::scale_color_manual(values = c("#1F78B4", "#33A02C")) +
    ggplot2::scale_x_continuous(
      labels = scales::unit_format(unit = "kb", scale = 1e-3),
      breaks = seq(-window.size, window.size, 1 * 10^(nchar(window.size - 1) - 1) * 2),
      limits = c(-window.size, window.size)
    ))
}
