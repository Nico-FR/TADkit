#' @title Plot 2 interaction matrices
#'
#' @description Same as `MATplot()` function but `mMATplot()` allow to plot 2 different matrices on the upper and lower part of the plot.
#' Two types of annotations can be added:
#' * domains (e.g. TADs or compartments): plot as triangles or lines on the upper or/and lower part of the matrix.
#' * interactions between domains/bins (e.g. loops): plot as rectangles on the upper and lower part of the matrix.
#'
#' @details The matrix input must be a `Matrix` or a `matrix` object for only one chromosome (see `cool2matrix()` function to read cool files).
#' All domains (TADs or compartments) are bed files (3 columns: chr, start and end) and can be R object (`dataframe` or `GRanges`) or the path of the files.
#' For `tad.lines`, another column can be used to specify different classes of domains (e.g compartment A or B). To use those domain classes, specify the column number (of the `tad.upper.line` and `tad.lower.line` inputs) with `tad.line.col` parameter and a custom set of colors with `line.colors` parameter.
#' Loop are stored in bedpe files (6 columns: chr1, start1, end, chr2, start2 and end2) and can be a `dataframe` object or the path of the file.
#' Chromosome domains and loops can be filter using `tad.chr` parameter.
#'
#' @inheritParams MATplot
#' @param matrix.upper,matrix.lower `dgCMatrix` or `matrix` object for only one chromosome.
#' @param matrix.upper.txt,matrix.lower.txt text to write on the upper or lower part of the matrix.
#'
#' @return `ggplot`
#'
#' @importFrom Matrix tril triu summary diag
#' @importFrom viridis scale_fill_viridis
#' @importFrom scales unit_format oob_squish_any
#' @importFrom dplyr full_join select
#' @importFrom BiocGenerics t
#' @import ggplot2
#' @export
#'
#' @examples
#' # get domains from boundaries:
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' domains.gr = boundary2domain(boundaries.gr)
#'
#' mMATplot(matrix.upper = mat_HCT116_chr19_50kb,
#'     matrix.lower = mat_HCT116_chr19_50kb,
#'     start = 5e6, stop = 15e6,
#'     bin.width = 50e3, log2 = TRUE,
#'     tad.upper.tri = domains.gr,
#'     tad.lower.tri = domains.gr,
#'     tad.chr = "chr19")
#'
#'  # add compartement A and B
#'  comp.gr = PC1calling(PC1_250kb.gr)
#'
#' mMATplot(matrix.upper = mat_HCT116_chr19_50kb,
#'     matrix.lower = mat_HCT116_chr19_50kb,
#'     start = 5e6, stop = 15e6,
#'     bin.width = 50e3, log2 = TRUE,
#'     tad.upper.tri = domains.gr,
#'     tad.lower.tri = domains.gr,
#'     tad.chr = "chr19",
#'     tad.upper.line = comp.gr,
#'     tad.lower.line = comp.gr,
#'     tad.line.col = 1)

mMATplot <- function(matrix.upper, matrix.lower, start, stop, bin.width, log2 = TRUE, scale.colors = "H", scale.limits = NULL,
                    matrix.upper.txt = NULL, matrix.lower.txt = NULL,
                    tad.upper.tri = NULL, tad.lower.tri = NULL, loop.bedpe = NULL, tad.chr = NULL, annotations.color = "red",
                    tad.upper.line = NULL, tad.lower.line = NULL, tad.line.col = NULL, line.colors = c("red", "blue")) {

  #local variables:
  i <- j <- x <- e <- s <- chr <- e2 <- s2 <-start1 <- end1 <- start2 <- end2 <- chr1 <- chr2 <- NULL

  #bin to read
  from = start %/% bin.width + 1 ; to = stop %/% bin.width #nb bin

  ####################
  #mat1
  if(!inherits(matrix.upper, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or Matrix object")}

  #filter matrix area
  mat1 = as(Matrix::triu(matrix.upper[from:to, from:to]), "CsparseMatrix")
  mat1[Matrix::triu(mat1 == 0)] <- NA

  #get log2
  if (log2 == T) {mat1@x = log2(mat1@x)}

  #squish
  if (!is.null(scale.limits)) {
    if (length(scale.limits) != 2) {stop("scale.limits must be a vector with 2 values")}
    mat1@x = scales::oob_squish_any(mat1@x, range = scale.limits)}

  #melt matrix
  upper_mat = Matrix::summary(Matrix::triu(mat1, 1))
  ####################

  ####################
  #mat2
  if(!inherits(matrix.lower, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or Matrix object")}

  #filter matrix area
  mat2 = as(Matrix::triu(matrix.lower[from:to, from:to]), "CsparseMatrix")
  mat2[Matrix::triu(mat2 == 0)] <- NA

  #get log2
  if (log2 == T) {mat2@x = log2(mat2@x)}

  #squish
  if (!is.null(scale.limits)) {
    if (length(scale.limits) != 2) {stop("scale.limits must be 2 value vector")}
    mat2@x = scales::oob_squish_any(mat2@x, range = scale.limits)}

  #melt matrix
  lower_mat = Matrix::summary(Matrix::tril(BiocGenerics::t(mat2), -1))

  ####################

  #merge matrices
  melted_mat = rbind(upper_mat, lower_mat)

  #add genomic coordinates
  melted_mat$j = (melted_mat$j + from - 1) * bin.width - bin.width / 2
  melted_mat$i = (melted_mat$i + from - 1) * - bin.width + bin.width / 2

  #geom_tile
  p = ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(y = i, x = j, fill = x))+
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
    ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())

  #scale_fill_gradient2
  if (scale.colors == "OE" | scale.colors == "ObsExp" | scale.colors == "OE2" | scale.colors == "ObsExp2") {
    p <- p + ggplot2::scale_fill_gradient2(
      low = ifelse(scale.colors %in% c("OE" ,"ObsExp"), "blue", "purple4"),
      high = ifelse(scale.colors %in% c("OE" ,"ObsExp"), "red", "darkgreen"),
      midpoint = 0, mid="white", na.value = "white")
  } else {
    p <- p + viridis::scale_fill_viridis(na.value = "black", option = scale.colors)
  }


  #upper tri
  if (!is.null(tad.upper.tri)) {

    if (is.character(tad.upper.tri)) {
      tad = read.table(tad.upper.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.upper.tri, "GRanges")) {
      tad = as.data.frame(tad.upper.tri, row.names = NULL)[,1:3]
    }

    if (is.data.frame(tad.upper.tri)) {
      tad = tad.upper.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = e2, yend = -s2), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = e2, y = -s2, xend = e, yend = -e), color = annotations.color, size = 0.3)
  }

  #lower tri
  if (!is.null(tad.lower.tri)) {

    if (is.character(tad.lower.tri)) {
      tad = read.table(tad.lower.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.lower.tri, "GRanges")) {
      tad = as.data.frame(tad.lower.tri, row.names = NULL)[,1:3]
    }

    if (is.data.frame(tad.lower.tri)) {
      tad = tad.lower.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -e2, xend = e, yend = -e), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = s2, yend = -e2), color = annotations.color, size = 0.3)
  }


  #upper line
  if (!is.null(tad.upper.line)) {

    if (is.character(tad.upper.line)) {
      tad = read.table(tad.upper.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.upper.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.upper.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.upper.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.upper.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
    }

    if (is.data.frame(tad.upper.line)) {
      tad = tad.upper.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    tad$col = as.factor(tad$col)
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)
      }

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -start, xend = e2, yend = -start, col = col), size = 1.5)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = stop, y = -s2, xend = stop, yend = -e2, col = col), size = 1.5)+
      ggplot2::scale_colour_manual(values = line.colors)
  }

  #lower line
  if (!is.null(tad.lower.line)) {

    if (is.character(tad.lower.line)) {
      tad = read.table(tad.lower.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.lower.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.lower.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.lower.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.lower.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
    }

    if (is.data.frame(tad.lower.line)) {
      tad = tad.lower.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    tad$col = as.factor(tad$col)
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)
      }

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = start, y = -s2, xend = start, yend = -e2, col = col), size = 1)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -stop, xend = e2, yend = -stop, col = col), size = 1)+
      ggplot2::scale_colour_manual(values = line.colors)
  }

  #loop
  if (!is.null(loop.bedpe)) {

    if (is.character(loop.bedpe)) {
      loop = read.table(loop.bedpe, header = F, sep = "\t")[,1:6]
    }

    if (is.data.frame(loop.bedpe)) {
      loop = loop.bedpe[,1:6]
    }

    names(loop) = c("chr1", "start1", "end1", "chr2", "start2", "end2")
    if (is.null(tad.chr)) {
      loop <- dplyr::filter(loop, start1 >= start, end1 < stop, start2 > start, end2 < stop)} else {
        loop <- dplyr::filter(loop, chr1 == tad.chr, chr2 == tad.chr, start1 >= start, end1 < stop, start2 > start, end2 < stop)}


    p = p + ggplot2::geom_rect(data = loop, ggplot2::aes(xmin = start2, xmax = end2, ymin = -start1, ymax = -end1), fill = NA, color = annotations.color, size = 0.3)+
      ggplot2::geom_rect(data = loop, ggplot2::aes(xmin = start1, xmax = end1, ymin = -start2, ymax = -end2), fill = NA, color = annotations.color, size = 0.3)
  }


  if (!is.null(matrix.upper.txt)) {
    p = p + annotate(geom = "text", x = stop, y = -start, label = matrix.upper.txt, col = "red", size = 10, hjust = "right", vjust = "top")
  }

  if (!is.null(matrix.lower.txt)) {
    p = p + annotate(geom = "text", x = start, y = -stop, label = matrix.lower.txt, col = "red", size = 10, hjust = "left", vjust = "left")
  }

  p
}





















