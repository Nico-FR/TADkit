#' @title Import matrix from HiC-Pro files
#'
#' @description `hicpro2matrix` imports the interaction matrix for one chromosome from HiC-Pro output files (.matrix and .bed, see details) as a `dgCMatrix` (upper triangular and sparse Matrix).
#' By default, it assumes the input file contains only the upper triangular part of the matrix.
#'
#' @details HiC-Pro is a tool for Hi-C data processing (https://github.com/nservant/HiC-Pro).
#' The `.matrix` file contains 3 columns: bin i, bin j, and the number of interactions.
#' The `.bed` file (index) contains 4 columns: chromosome, start, end, and bin ID.
#' `hicpro2matrix` filters the interactions to keep only intra-chromosomal interactions for the specified chromosome.
#' The resulting matrix is upper triangular.
#'
#' @param matrix.path The path to the HiC-Pro .matrix file.
#' @param index.path The path to the HiC-Pro .bed index file.
#' @param chr The selected chromosome.
#' @param verbose Logical. Whether or not to print messages. Default = `TRUE`.
#'
#' @return A `dgCMatrix` object: upper triangular and sparse Matrix.
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table fread :=
#' @export
#'
hicpro2matrix <- function(matrix.path, index.path, chr, verbose = TRUE) {

  # Use a different name for the variable to avoid conflict with data.table column names
  selected_chr <- chr

  # 1. Read index file
  if (verbose) message("Reading index file...")
  # HiC-Pro index: chr, start, end, id
  index <- data.table::fread(index.path, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "id"))

  # 2. Filter index for the chromosome
  index_chr <- index[index[["chr"]] == selected_chr, ]
  if (nrow(index_chr) == 0) {
    stop(paste0("Chromosome '", selected_chr, "' not found in index file."))
  }

  min_id <- min(index_chr$id)
  max_id <- max(index_chr$id)
  nb_bins <- nrow(index_chr)

  if (verbose) message(paste0("Parsing interactions for ", selected_chr, " (", nb_bins, " bins)..."))

  # 3. Read matrix file
  # The file has 3 columns: bin1_id, bin2_id, count
  # We use colClasses for speed and memory efficiency.
  if (verbose) message("Reading matrix file...")
  mat_data <- data.table::fread(matrix.path,
                                header = FALSE,
                                col.names = c("bin1_id", "bin2_id", "count"),
                                colClasses = c("integer", "integer", "numeric"))

  # 4. Filter for intra-chromosomal interactions
  # Keep only if both bins are within the chromosome's range
  if (verbose) message("Filtering interactions for the selected chromosome...")
  mat_data <- mat_data[bin1_id >= min_id & bin1_id <= max_id &
                       bin2_id >= min_id & bin2_id <= max_id]

  if (nrow(mat_data) == 0) {
    stop("No interactions found for this chromosome.")
  }

  # 5. Re-index and keep only upper triangular part
  # Re-index to start from 1 for the specific chromosome
  mat_data[, `:=`(bin1_id = bin1_id - min_id + 1, bin2_id = bin2_id - min_id + 1)]

  # Filter to ensure it is upper triangular (as per assumption)
  mat_data <- mat_data[bin1_id <= bin2_id, list(bin1_id, bin2_id, count)]

  # 6. Create sparse matrix
  m <- Matrix::sparseMatrix(i = mat_data$bin1_id,
                            j = mat_data$bin2_id,
                            x = mat_data$count,
                            dims = c(nb_bins, nb_bins))

  return(m)
}
