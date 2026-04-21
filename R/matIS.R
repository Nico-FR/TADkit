#' Calculate Insulation Score (IS) on matrices
#'
#' This function calculates the Insulation Score stabilized with log1p.
#'
#' @param matrix.lst List of `dgCMatrix` or `matrix` objects for a single chromosome.
#' @param bin.width Integer. The matrix bin width in base pairs (e.g., 100000).
#' @param square_size Integer. The size of the sliding square window (default: 300kb).
#' @param start Integer. Start coordinate of the region of interest (optional).
#' @param stop Integer. Stop coordinate of the region of interest (optional).
#' @param seqname Character. Chromosome name, default = "NA".
#'
#' @details
#' The Insulation Score (IS) is calculated for each bin as the mean of the interaction counts in a square window centered on the bin.
#' The final IS values are normalized the mean IS across the analyzed region. Therefore considering the entire chromosome for better normalization is recommended.
#'
#' @examples
#' matIS(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb),
#'   bin.width = 50e3, square_size = 300e3,
#'   start = 5e6, stop = 15e6, seqname = "chr19")
#'
#' @return A list of `GRanges` objects with IS scores.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
matIS <- function(matrix.lst, bin.width, square_size = 300e3, start = NULL, stop = NULL, seqname = "NA") {

  ################################
  # Sanity checks
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  n_rows <- sapply(matrix.lst, nrow)
  n_cols <- sapply(matrix.lst, ncol)
  if (length(unique(n_rows)) != 1 || length(unique(n_cols)) != 1) {
    stop("All matrices in matrix.lst must have the same sizes")
  }

  if (is.null(names(matrix.lst))) {
    names(matrix.lst) <- 1:length(matrix.lst)
  }
  ################################

  ################################
  # Matrix parameters
  squareBin <- round(square_size / bin.width)
  if (round(square_size / bin.width) != square_size / bin.width) {
    message("square_size is not a multiple of bin.width, rounded to ", squareBin * bin.width, " (i.e. ", squareBin, " bins).")
  }
  n_bins_total <- unique(n_rows)[1]

  # Matrix coordinate of limits
  l.start <- 1
  l.stop  <- n_bins_total

  if (!is.null(start)) {
    tmp <- floor(start / bin.width + 1) - squareBin
    l.start <- max(1, tmp)
  }
  if (!is.null(stop)) {
    tmp <- ceiling(stop / bin.width) + squareBin
    l.stop <- min(n_bins_total, tmp)
  }

  n_bins_roi <- l.stop - l.start + 1 # Number of bins in the cropped area
  ################################

  array_IS.lst <- list()

  for (mat_index in 1:length(matrix.lst)) {

    mat <- matrix.lst[[mat_index]][l.start:l.stop, l.start:l.stop]
    array_IS <- rep(NA, n_bins_roi)

    # Compute IS for each bin
    for (i in 1:n_bins_roi) {

      # 1. Edge filtering: skip if the square exceeds matrix boundaries
      if ((i - squareBin < 1) || (i + squareBin > n_bins_roi)) {
        next
      }

      # 2. Define the square limits
      # Upstream limits (rows)
      start_U <- i - squareBin
      end_U   <- i - 1

      # Downstream limits (columns)
      start_D <- i + 1
      end_D   <- i + squareBin

      # 3. Extract the square submatrix and compute the mean
      square_mat <- array(mat[start_U:end_U, start_D:end_D])
      score <- mean(square_mat, na.rm = TRUE)

      # 4. Strict safety filters (matching the Python script)
      if (is.na(score) || score == 0) {
        next
      }

      array_IS[i] <- score
    }

    # 5. Normalization step (log1p of the ratio to the global mean)
    global_mean <- mean(array_IS, na.rm = TRUE)
    array_IS <- log1p(array_IS / global_mean)

    array_IS.lst[[mat_index]] <- array_IS
  }

  ################################
  # Merge output list into GRanges
  start1 <- ifelse(is.null(start), (l.start - 1) * bin.width, start)
  stop1  <- ifelse(is.null(stop), l.stop * bin.width, stop)

  df.tmp <- data.frame(
    seqname = rep(seqname, n_bins_roi),
    start   = ((l.start:l.stop) - 1) * bin.width,
    stop    = (l.start:l.stop) * bin.width
  )

  bedgraph.lst <- lapply(1:length(array_IS.lst), function(i) {
    gr_df <- cbind(df.tmp, IS = array_IS.lst[[i]])
    gr_df <- gr_df[gr_df$stop > start1 & gr_df$start < stop1, ]
    GenomicRanges::makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  })

  names(bedgraph.lst) <- names(matrix.lst)
  return(bedgraph.lst)
}
