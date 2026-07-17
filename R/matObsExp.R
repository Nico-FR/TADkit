#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction count observed) `matObsExp()` return the ratio: observed / expected.
#' The input matrix must be a `Matrix` or `matrix` object (i.e matrix with as many rows and columns than the number of bin).
#' The output can be plot with `MATplot(output, log2 = TRUE, scale.colors = "ObsExp")`.
#'
#' @details The expected number of interaction corresponds to the average interaction counts according to bin distances (only non-zero values are used to compute).
#' Note that the expected number of interaction is only estimated from the chromosome supplied. Other genome wide approaches may be considered.
#'
#' @param matrix `Matrix` or `matrix` object.
#' @param output default is "OE" to return observed / expected matrix. Use "E" or "Exp" to return expected matrix.
#'
#' @return `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @importFrom stats toeplitz
#' @importFrom Matrix triu summary sparseMatrix
#' @importFrom dplyr mutate filter group_by summarise
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#' @importFrom tidyr complete
#'
#' @export
#'
#' @examples
#' mat_obsexp = matObsExp(mat_HCT116_chr19_50kb)
#' MATplot(matrix = mat_obsexp,
#'     start = 5e6, stop = 15e6,
#'     bin.width = 50e3,
#'     log2 = TRUE, scale.colors = "OE")
#'

matObsExp <- function(matrix, output = "OE") {

  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  if(inherits(matrix, "matrix")) {
    matrix = methods::as(matrix, "CsparseMatrix")}

  n = ncol(matrix)

  s_original = Matrix::summary(matrix)
  s = s_original[s_original$i <= s_original$j, ]
  s$dist = s$j - s$i

  diag_mean_vector = rep(NaN, n)
  if (nrow(s) > 0) {
    means = tapply(s$x, s$dist, mean, na.rm = TRUE)
    dists = as.integer(names(means))
    diag_mean_vector[dists + 1] = means
  }

  if (output == "OE") {
    s_full = s_original
    s_full$dist = abs(s_full$i - s_full$j)
    s_full$x = s_full$x / diag_mean_vector[s_full$dist + 1]

    out = Matrix::sparseMatrix(
      i = s_full$i,
      j = s_full$j,
      x = s_full$x,
      dims = dim(matrix),
      symmetric = FALSE
    )
  } else {
    mat_expected = stats::toeplitz(diag_mean_vector)
    out = mat_expected
  }

  return(
    if(inherits(out, "CsparseMatrix")) {out} else {methods::as(out, "CsparseMatrix")}
  )
}









