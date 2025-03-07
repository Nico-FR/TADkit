#' @title Pileup of matrices surrounding genomic features
#'
#' @description This function allows the quantification of patterns (of the matrix) surrounding genomic features.
#' For example this function allow to visualize the patterns of the matrices around TADs boundary.
#'
#' @details From genomic features (e.g. start of TADs):
#'  - export the matrix +/- `window.size` of each genomic features,
#'  - pileup the matrices (i.e. sum all the matrices),
#'  - return the pileup matrix or plot the observed/expected of the pileup matrix.
#'
#' @inheritParams MATplot
#' @param annot.gr `GRanges` with genomic annotations.
#' @param chr The selected chromosome used to filter `annot.gr`.
#' @param annot.boundary Type of feature to analyzed. `"Start"`, `"end"` or `"center"` of each `annot.gr`.
#' @param window.size Window to analyze the matrix on each side of the features. By default the window is 40 times the `bin.width`.
#' @param output Default is `"matrix"` to return the `matrix` (i.e sum of matrices). Use `"plot"` to return a `ggplot` of the observed / expected of the pileup of the matrices, i.e. ObsExp(sum of matrices).

#'
#' @return A `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @examples
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' domains.gr = boundary2domain(boundaries.gr)
#' MATfeatures(matrix = mat_HCT116_chr19_50kb, bin.width = 50e3, annot.gr = domains.gr, chr = "chr19",
#'             annot.boundary = "start", window.size = 1e6, output = "plot")
#'
#' @export
#'
MATfeatures <- function(matrix, bin.width, annot.gr, chr, annot.boundary = "start", window.size = NULL, output = "matrix") {

  i <- j <- x <- NULL

  ############################################################
  #sanity check
  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or Matrix object")}

  if (is.null(window.size)) {window.size = 40 * bin.width}
  ############################################################

  #features positions
  features.gr <- GenomicRanges::resize(annot.gr, 1, fix = annot.boundary) %>%
    #removed features too close (ie distance < window.size) to the border of the matrix
    GenomicRanges::restrict(start = window.size + bin.width + 1, end = nrow(matrix) * bin.width - window.size -bin.width - 1)

  #bin ranges of the matrix
  bin.matrix.gr = GenomicRanges::GRanges(seqnames = chr,
                                         ranges = IRanges::IRanges(start = 0:(nrow(matrix) - 1) * bin.width + 1,
                                                                   end = 1:(nrow(matrix)) * bin.width))


  bin.matched.lst = S4Vectors::subjectHits(GenomicRanges::findOverlaps(features.gr, bin.matrix.gr))

  message(paste0("Staking of ", length(bin.matched.lst), " matrices on chr ", chr, "."))

  #if no overlap: return empty matrix else: return stacked matrix around each feature
  if (identical(bin.matched.lst, integer(0))) {

    warning(paste0("no features on chromosome ", chr))

    pil_mat = matrix(0, nrow = (window.size %/% bin.width * 2) + 2, ncol = (window.size %/% bin.width * 2) + 2)
  } else {

    mat.lst = sapply(bin.matched.lst, function(x){

      #bin to read
      from = x - window.size %/% bin.width #upstream bin
      to = x + window.size %/% bin.width + 1 #downstream bin

      #filter matrix area
      Matrix::triu(matrix[from:to, from:to])
    })

    #staking matrices
    pil_mat = base::Reduce(`+`, mat.lst) %>% as.matrix %>% methods::as("CsparseMatrix")
  }

  if (output == "matrix") {return(pil_mat)}

  if (output == "plot") {

    mat = TADkit::matObsExp(pil_mat)
    mat[Matrix::triu(mat == 0)] <- NA
    mat@x = log2(mat@x)

    upper_mat = Matrix::summary(Matrix::triu(mat, 1))
    diag_mat = data.frame(i = 1:nrow(mat), j = 1:nrow(mat), x = Matrix::diag(mat))
    lower_mat = data.frame(i = upper_mat$j, j = upper_mat$i, x = upper_mat$x)

    melted_mat = rbind(upper_mat, lower_mat, diag_mat)
    melted_mat$j = (melted_mat$j - 0.5) * bin.width - bin.width / 2 - window.size
    melted_mat$i = (melted_mat$i - 0.5) * -bin.width + bin.width / 2 + window.size
    p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(y = i, x = j, fill = x))+
      ggplot2::scale_fill_gradient2(low = "blue", high ="red",midpoint = 0, mid="white", na.value = "white")+
      ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
      ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())+
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)+geom_vline(xintercept = 0, linetype = "dashed", size = 0.25)

    return(p)
  }
}
