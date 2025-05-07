#' @title Quantify interactions changes for each bin between 2 matrices
#'
#' @description For each bin, two metrics can be computed to quantify interactions changes between 2 matrices (mutated and wildtype matrices) :
#'
#' - the correlation of interaction values,
#'
#' - the absolute mean fold change of interaction values.
#'
#' @param mutated.mat,wildtype.mat mutated (query) or wildtype matrices as `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Bin width of the matrix in base pair.
#' @param start,stop Area in bp of the chromosome to be computed. Default is NULL to use the entire chromosome (i.e. entire matrix).
#' @param max.distance Maximum distance between interactions. Default is 2e6 (i.e. 2Mb): to keep only interactions between 2 bins with a distance smaller than 2Mb.
#' @param method Default is "Corr" for pearson correlation. Use "FC" to compute the mean fold change between interactions.
#' @param seqname chromosome names as character for output, default = "1".
#' @param log2 Default is TRUE, to use log2 of the interaction values before computing correlations. Necessary if interaction values are not normalized by distance (Observed over Expected).
#'
#' @importFrom magrittr %>%
#'
#' @return GRanges
#'
#' @export
#'
#' @examples
#' # to do
#'
binDiff <- function(mutated.mat, wildtype.mat, bin.width, start = NULL, stop = NULL, max.distance = 2e6, method = "Corr", seqname = "1", log2 = TRUE) {

   #mutated.mat = mat.lst[[1]] ; wildtype.mat = mat.lst[[2]] ; bin.width = 64e3
   #start = NULL; stop = NULL;  method = "Corr"; seqname = "4"; log2 = TRUE
   #max.distance = 640e3
  #start = (65642800 %/% bin.width + 1) * bin.width
  #stop = (76557789 %/% bin.width) * bin.width

  if (is.null(bin.width)) {
    stop("bin.width value is needed!")
  }

  # Convert input matrices to list of matrices
  matrix.lst = list(as.matrix(mutated.mat), as.matrix(wildtype.mat))

  # Filter area based on start and stop parameters
  if(!(is.null(start) | is.null(stop))) { #if start & stop != NULL

    # Calculate bin indexes
    if (start %% bin.width == 0 & stop %% bin.width == 0) {

      if (start > stop) {
        stop("start cannot be greater than stop!")
      }

      start = start / bin.width + 1
      stop = stop / bin.width
      } else {
        start = round(start / bin.width) + 1 # +1 because matrix coordinates start to 1. So first nucleotide is bin nb 1.
        stop = ifelse((round(stop / bin.width) + 1) <= start, start, round(stop / bin.width))
        message(paste0("start/stop are not multiples of bin.width, round to ", (start - 1) * bin.width,
                       " and ", stop * bin.width, " (i.e. ", stop + 1 - start, " bins)."))
      }

    # Filter matrices based on start and stop indexes
    matrix.lst = lapply(matrix.lst, function(mat){
      mat[start:stop,start:stop]
    })
  } else {
    start = 1
    stop = nrow(matrix.lst[[1]])
  }

  # Ensure matrices are symmetric and replace zeros with NA
  matrix.lst = lapply(matrix.lst, function(mat) {
    if (!isSymmetric(mat)) {
      mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    }
    mat[mat == 0] <- NA
    return(mat %>% as.matrix)
  })

  # Nb bin in the matrices
  nb_row = nrow(matrix.lst[[1]])

  # Compute nb bin according to max.distance
  if(!is.null(max.distance)) {
    bin_max_dist = max.distance %/% bin.width
  } else {
    bin_max_dist = nb_row
  }

  #compute score for each bin
  score <- sapply(1:nb_row, function(i) {

    #compute bin indexes to keep
    start_i <- max(1, i - bin_max_dist + 1)
    stop_i <- min(nb_row, i + bin_max_dist - 1)

    # Get the i-th row of the first matrix and filter interactions by bin_max_dist (i.e. max.distance)
    row1 <- matrix.lst[[1]][i,start_i:stop_i]
    # Get the i-th row of the second matrix
    row2 <- matrix.lst[[2]][i,start_i:stop_i]

    if (method == "Corr") {
      # Apply log2 transformation if log2 is TRUE
      if(isTRUE(log2)) {
        row1 <- log2(row1)
        row2 <- log2(row2)
      }

      # Calculate the correlation between the two rows
      stats::cor(row1, row2, use = "pairwise.complete.obs")
    } else if (method == "FC") {
      # Calculate the fold change between the two rows
      (row1 / row2) %>% log2 %>% abs %>% mean(na.rm = TRUE)
      }
    })

  gr = data.frame(
    seqnames = seqname,
    start = ((start:stop) - 1) * bin.width,
    stop = ((start:stop)) * bin.width,
    score = score
    ) %>%
    `colnames<-`(c("seqnames", "start", "end", method)) %>% # Set column names
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) # Convert to GRanges object

  return(gr)
}




