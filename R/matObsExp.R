#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction count observed) `matObsExp()` return the ratio: observed / expected.
#' The input matrix must be a `Matrix` or `matrix` object (i.e matrix with as many rows and columns than the number of bin).
#' The output can be plot with `MATplot(output, log2 = TRUE, scale.colors = "ObsExp")`.
#'
#' @details The expected number of interaction corresponds to the average interaction counts according to bin distances.
#' Note that the expected number of interaction is only estimated from the chromosome supplied. Other genome wide approaches may be considered.
#'
#' @param matrix `Matrix` or `matrix` object.
#' @param output default is "OE" to return observed / expected matrix. Use "E" or "Exp" to return expected matrix.
#'
#' @return `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @importFrom stats toeplitz
#' @importFrom Matrix triu summary
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

  i <- j <- . <- dist <- x <- NULL

  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  if(inherits(matrix, "matrix")) {
    matrix = methods::as(matrix, "CsparseMatrix")}

  diag_mean.df = Matrix::summary(matrix) %>%
    dplyr::mutate(i = factor(i, levels = 1:ncol(matrix))) %>% #add missing index as levels
    dplyr::mutate(j = factor(j, levels = 1:ncol(matrix))) %>% #add missing index as levels
    tidyr::complete(., i, j, fill = list(x = NA), explicit = FALSE) %>% #add value NA to missing bins
    dplyr::filter(as.numeric(i) <= as.numeric(j)) %>% #filter lower matrix
    dplyr::mutate(dist = abs(as.numeric(i) - as.numeric(j))) %>% #add distances between bins
    dplyr::group_by(dist) %>% dplyr::summarise(diag_mean = mean(x, na.rm = TRUE)) # mean according to distances
  mat_expected = stats::toeplitz(diag_mean.df$diag_mean)

  out = if(output == "OE") {matrix / mat_expected} else if(output == "E") {mat_expected} else if(output == "Exp") {mat_expected}

  return(
    if(inherits(out, "CsparseMatrix")) {out} else {methods::as(out, "CsparseMatrix")}
  )
}












