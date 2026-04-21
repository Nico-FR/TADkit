#' Calculate Directionality Index (DI) on matrices
#'
#' This function calculates the DI (stabilized with log1p).
#'
#' @param matrix.lst List of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Integer. The matrix bin.width in base pairs (e.g., 100000).
#' @param DI_distance Integer. The observation window distance (default is 500kb). It is recommended to use a distance between 400kb to 1mb according to bin.width (e.g. 500kb for bin.width = 16kb, 1mb for bin.width = 128kb).
#' @param start Integer. Start coordinate of the region of interest in bp (optional).
#' @param stop Integer. Stop coordinate of the region of interest in bp (optional).
#' @param seqname Character. Chromosome name, default = "NA".
#'
#' @return A list of `GRanges` objects with DI scores.
#' @export
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @examples
#' matDI(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb),
#'   bin.width = 50e3, DI_distance = 500e3,
#'   start = 5e6, stop = 15e6, seqname = "chr19")
#'
matDI <- function(matrix.lst, bin.width, DI_distance = 500e3, start = NULL, stop = NULL, seqname = "NA") {

  ################################
  # Sanity checks
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  # Check if all matrices have the same size
  n_rows <- sapply(matrix.lst, nrow)
  n_cols <- sapply(matrix.lst, ncol)
  if (length(unique(n_rows)) != 1 || length(unique(n_cols)) != 1) {
    stop("All matrices in matrix.lst must have the same sizes")
  }

  # Add names to matrix.lst if not already named
  if (is.null(names(matrix.lst))) {
    names(matrix.lst) <- 1:length(matrix.lst)
  }
  ################################

  ################################
  # Matrix parameters
  distanceBin <- round(DI_distance / bin.width) # DI_distance in bins
  if (round(DI_distance / bin.width) != DI_distance / bin.width) {
    message("DI_distance is not a multiple of bin.width, rounded to ", distanceBin * bin.width, " (i.e. ", distanceBin, " bins).")
  }
  n_bins_total <- unique(n_rows)[1] # Total number of bins in the matrix

  # Matrix coordinate of limits
  l.start <- 1
  l.stop  <- n_bins_total

  if (!is.null(start)) {
    tmp <- floor(start / bin.width + 1) - distanceBin
    l.start <- max(1, tmp)
  }
  if (!is.null(stop)) {
    tmp <- ceiling(stop / bin.width) + distanceBin
    l.stop <- min(n_bins_total, tmp)
  }

  n_bins_roi <- l.stop - l.start + 1 # Number of bins in the cropped area
  ################################

  array_DI.lst <- list() # Initialize the output list

  for (mat_index in 1:length(matrix.lst)) {

    mat <- matrix.lst[[mat_index]][l.start:l.stop, l.start:l.stop] # Crop matrix
    array_DI <- rep(NA, n_bins_roi) # Initialize the output vector with NA

    # Compute DI for each bin
    for (i in 1:n_bins_roi) {

      # 1. Edge filtering
      if ((i - distanceBin < 1) || (i + distanceBin > n_bins_roi)) {
        next
      }

      # 2. Calculate sum U (Upstream / Left interactions)
      start_U <- i - distanceBin
      end_U   <- i - 1
      U <- sum(mat[start_U:end_U, i], na.rm = TRUE)

      # 3. Calculate sum D (Downstream / Right interactions)
      start_D <- i + 1
      end_D   <- i + distanceBin
      D <- sum(mat[i, start_D:end_D], na.rm = TRUE)

      # Expected value (mean of U and D)
      E <- (U + D) / 2

      # 4. Strict safety filters
      if (is.na(U + D) || (D - U == 0) || min(U, D) == 0) {
        next
      }

      # 5. Directionality Sign
      sign_val <- (D - U) / abs(D - U)

      # 6. Chi-square component with Log1p stabilization
      chi_sq <- ((U - E)^2) / E + ((D - E)^2) / E
      array_DI[i] <- sign_val * log1p(chi_sq)
    }

    array_DI.lst[[mat_index]] <- array_DI
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

  bedgraph.lst <- lapply(1:length(array_DI.lst), function(i) {
    gr_df <- cbind(df.tmp, DI = array_DI.lst[[i]])
    gr_df <- gr_df[gr_df$stop > start1 & gr_df$start < stop1, ]
    GenomicRanges::makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
  })

  names(bedgraph.lst) <- names(matrix.lst)
  return(bedgraph.lst)
}
