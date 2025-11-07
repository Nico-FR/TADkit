#' @title Plot intra-chromosomal interaction matrix
#'
#' @description `interMATplot` allow to plot non-symmetrical matrix with 2 types of annotations:
#' * domains (e.g. TADs or compartments): plot as triangles and/or lines on the x or/and y axis of the matrix.
#' * interactions between domains/bins (e.g. loops): plot as squares or rectangles.
#'
#' @details The matrix input must be a `Matrix` or a `matrix` object with interactions from 2 chromosome (i.e. inter-chromosomal interaction).
#' All domains (TADs or compartments) are bed files (3 columns: chr, start and end) and can be R object (`dataframe` or `GRanges`) or the path of the files.
#' For `tad.lines`, another column can be used to specify different classes of domains (e.g compartment A or B). To use those domain classes, specify the column number (from `tad.x.line` and `tad.y.line` inputs) with `tad.line.col` parameter and a custom set of colors with `line.colors` parameter.
#' Loop are stored in bedpe files (6 columns: chr_y, start_y, end_y, chr_x, start_x and end_x) and can be a `dataframe` object or the path of the file.
#' Chromosome domains and loops can be filter using `tad.chr.y` and `tad.chr.y` parameter.
#'
#' @param matrix non-symmetrical `dgCMatrix` or `matrix` object with inter-chromosomal interactions.
#' @param bin.width Bin width of the matrix in base pair.
#' @param start_y,stop_y Region of interest in base pair on y-axis. Default is NULL to plot all the chromosome.
#' @param start_x,stop_x Region of interest in base pair on x-axis. Default is NULL to plot all the chromosome.
#' @param log2 logical. Use the log2 of the matrix values. Default is `TRUE`.
#' @param scale.limits Use to set limits on the scale. Default is NULL to use all values. E.g if use scale.limits = c(0, 10): Values < 0 will be fix to 0 and values > 10 will be fix to 10.
#' @param scale.colors A character string indicating the color map option to use. Eight colors palettes are available from `viridis` package. Two supplementary palettes `"OE"` and  `"OE2"` (blue to red and purple to green respectively) are made for data centered on 0 (e.g log2 of observed/expected matrix). Default is `"H"`:
#' * `"magma"` (or `"A"`),
#' * `"inferno"` (or `"B"`),
#' * `"plasma"` (or `"C"`),
#' * `"viridis"` (or `"D"`),
#' * `"cividis"` (or `"E"`),
#' * `"rocket"` (or `"F"`),
#' * `"mako"` (or `"G"`),
#' * `"turbo"` (or `"H"`),
#' * `"ObsExp"` (or `"OE"`),
#' * `"ObsExp2"` (or `"OE2"`).
#' @param tad.y.tri,tad.x.tri `data.frame`, `GRanges` or the bed files path with the TAD to plot as triangle in the lower (x) or upper (y) part of the diagonal. Default is `NULL`.
#' @param tad.y.line,tad.x.line `data.frame`, `GRanges` or the bed files path with the TAD to plot as line in the x or y axis of the matrix. Default is `NULL`.
#' @param tad.line.col Column number of `tad.y.line` and `tad.x.line` files that contain factors used to color lines. Default is `NULL`.
#' @param loop.bedpe `data.frame` or bedpe files path to plot on the matrix. Six columns table (chr_y, start_y, end_y, chr_x, start_x, end_x) that gives areas between 2 chromosomes. Default is `NULL`
#' @param tad.chr.y,tad.chr.x Chromosome name to filter annotations (domains and loop). Default is `NULL`
#' @param annotations.color Color for loop and triangular annotations. Default is `"red"`.
#' @param line.colors Colors for `tad.x.line` and `tad.y.line`. Default is `c("red", "blue")`.
#'
#' @return `ggplot`
#'
#' @importFrom Matrix triu summary diag tril
#' @importFrom viridis scale_fill_viridis
#' @importFrom scales unit_format oob_squish_any
#' @importFrom dplyr filter
#' @importFrom BiocGenerics t
#' @import ggplot2
#' @export
#'
#' @examples
#'
#' # create a matrix 500x1000 bins
#' matrix <- matrix(rnorm(500000, mean = 50, sd = 10) + c(rep(0, 250), seq(0, 50, length.out = 500), rep(0, 250)), 500, 1000)
#'
#' interMATplot(matrix,
#'     bin.width = 100e3)
#'
interMATplot <- function(matrix, bin.width, start_y = NULL, stop_y = NULL, start_x = NULL, stop_x = NULL, log2 = T, scale.colors = "H", scale.limits = NULL,
                    tad.x.tri = NULL, tad.y.tri = NULL, loop.bedpe = NULL, tad.chr_x = NULL, tad.chr_y = NULL, annotations.color = "red",
                    tad.x.line = NULL, tad.y.line = NULL, tad.line.col = NULL, line.colors = c("red", "blue")) {

  #matrix = cool2matrix(cattle.files$norm_mcool, chr = 25, chr2 = 26, bin.width = 128000,
   #                    balance = FALSE, verbose = TRUE)
  #bin.width = 128e3; start_y = NULL; stop_y = 40e6; start_x = NULL; stop_x = 50e6; log2 = T; scale.colors = "H"; scale.limits = NULL;
  #tad.x.tri = NULL; tad.y.tri = NULL; loop.bedpe = NULL; tad.chr_x = NULL; tad.chr_y = NULL; annotations.color = "red";
  #tad.x.line = NULL; tad.y.line = NULL; tad.line.col = NULL; line.colors = c("red", "blue")
  #tad.x.line = cattle.files$dchic_comp; tad.chr_x = 26
  #tad.y.line = cattle.files$dchic_comp; tad.chr_y = 25

  #sanity check
  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  #starts and stops if NULLs
  if (is.null(start_y)) {start_y = 1}
  if (is.null(stop_y)) {stop_y = nrow(matrix) * bin.width}
  if (is.null(start_x)) {start_x = 1}
  if (is.null(stop_x)) {stop_x = ncol(matrix) * bin.width}

  # starts and stops in bin indexes
  from1 = start_y %/% bin.width + 1 ; to1 = stop_y %/% bin.width
  from2 = start_x %/% bin.width + 1 ; to2 = stop_x %/% bin.width

  #filter matrix area
  mat = as(Matrix::Matrix(matrix[from1:to1, from2:to2]), "CsparseMatrix")
  mat[Matrix::Matrix(mat == 0)] <- NA # set 0 to NA

  #get log2
  if (log2 == T) {mat@x = log2(mat@x)}

  #squish
  if (!is.null(scale.limits)) {
    if (length(scale.limits) != 2) {stop("scale.limits must be a vector with 2 values")}
    mat@x = scales::oob_squish_any(mat@x, range = scale.limits)}

  #melt matrix
  melted_mat = Matrix::summary(mat) %>%
    dplyr::rename(c("y" = "i", "x" = "j", "i" = "x")) %>% #rename to avoid confusion
    dplyr::mutate(x = (x + from2 - 1.5) * bin.width) %>%  # genomic coordinates = center of the bin (needed for geom_tile)
    dplyr::mutate(y = -(y + from1 - 1.5) * bin.width)

  #geom_tile
  p = ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(y = y, x = x, fill = i))+
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6),
                                limits = c(min(melted_mat$x - bin.width / 2), max(melted_mat$x + bin.width / 2)))+ #limit to match geom_tile centers
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6),
                                limits = c(max(melted_mat$y + bin.width / 2), min(melted_mat$y - bin.width / 2)))+
    ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())

  #scales colors
  if (scale.colors == "OE" | scale.colors == "ObsExp" | scale.colors == "OE2" | scale.colors == "ObsExp2") {
    p <- p + ggplot2::scale_fill_gradient2(
      low = ifelse(scale.colors %in% c("OE" ,"ObsExp"), "blue", "purple4"),
      high = ifelse(scale.colors %in% c("OE" ,"ObsExp"), "red", "darkgreen"),
      midpoint = 0, mid="white", na.value = "white")
  } else {
    p <- p + viridis::scale_fill_viridis(na.value = "black", option = scale.colors)
    }

  #tad.x.tri (horizontal axes)
  if (!is.null(tad.x.tri)) {

    if (is.character(tad.x.tri)) {
      tad = read.table(tad.x.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.x.tri, "GRanges")) {
      tad = as.data.frame(tad.x.tri, row.names = NULL)[,1:3]
    }

    if (is.data.frame(tad.x.tri)) {
      tad = tad.x.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr_x)) {
      tad <- dplyr::filter(tad, e > start_x, s < stop_x)} else {
        tad <- dplyr::filter(tad, chr == tad.chr_x, e > start_x, s < stop_x)}

    tad$e2 <- ifelse(tad$e >= stop_x, stop_x, tad$e)
    tad$s2 <- ifelse(tad$s <= start_x, start_x, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = e2, yend = -s2), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = e2, y = -s2, xend = e, yend = -e), color = annotations.color, size = 0.3)
  }

  #tad.y.tri (vertical axes)
  if (!is.null(tad.y.tri)) {

    if (is.character(tad.y.tri)) {
      tad = read.table(tad.y.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.y.tri, "GRanges")) {
      tad = as.data.frame(tad.y.tri, row.names = NULL)[,1:3]
    }

    if (is.data.frame(tad.y.tri)) {
      tad = tad.y.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr_y)) {
      tad <- dplyr::filter(tad, e > start_y, s < stop_y)} else {
        tad <- dplyr::filter(tad, chr == tad.chr_y, e > start_y, s < stop_y)}

    tad$e2 <- ifelse(tad$e >= stop_y, stop_y, tad$e)
    tad$s2 <- ifelse(tad$s <= start_y, start_y, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -e2, xend = e, yend = -e), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = s2, yend = -e2), color = annotations.color, size = 0.3)
  }


  #tad.x.line (horizontal lines)
  if (!is.null(tad.x.line)) {

    if (is.character(tad.x.line)) {
      tad = read.table(tad.x.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.x.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.x.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.x.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.x.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
    }

    if (is.data.frame(tad.x.line)) {
      tad = tad.x.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    tad$col = as.factor(tad$col)
    if (is.null(tad.chr_x)) {
      tad <- dplyr::filter(tad, e > start_x, s < stop_x)} else {
        tad <- dplyr::filter(tad, chr == tad.chr_x, e > start_x, s < stop_x)
      }

    tad$e2 <- ifelse(tad$e >= stop_x, stop_x, tad$e)
    tad$s2 <- ifelse(tad$s <= start_x, start_x, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -start_y, xend = e2, yend = -start_y, col = col), size = 1.5)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -stop_y, xend = e2, yend = -stop_y, col = col), size = 1.5)
  }

  #tad.y.line
  if (!is.null(tad.y.line)) {

    if (is.character(tad.y.line)) {
      tad = read.table(tad.y.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.y.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.y.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.y.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.y.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
    }

    if (is.data.frame(tad.y.line)) {
      tad = tad.y.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    tad$col = as.factor(tad$col)
    if (is.null(tad.chr_y)) {
      tad <- dplyr::filter(tad, e > start_y, s < stop_y)} else {
        tad <- dplyr::filter(tad, chr == tad.chr_y, e > start_y, s < stop_y)
      }

    tad$e2 <- ifelse(tad$e >= stop_y, stop_y, tad$e)
    tad$s2 <- ifelse(tad$s <= start_y, start_y, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = start_x, y = -s2, xend = start_x, yend = -e2, col = col), size = 1)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = stop_x, y = -s2, xend = stop_x, yend = -e2, col = col), size = 1)
  }

  if (!is.null(tad.x.line) | !is.null(tad.y.line)) {
    p = p + ggplot2::scale_color_manual(values = line.colors)
  }

  #loop
  if (!is.null(loop.bedpe)) {

    if (is.character(loop.bedpe)) {
      loop = read.table(loop.bedpe, header = F, sep = "\t")[,1:6]
    }

    if (is.data.frame(loop.bedpe)) {
      loop = loop.bedpe[,1:6]
    }

    names(loop) = c("c_y", "s_y", "e_y", "c_x", "s_x", "e_x")
    if (is.null(tad.chr)) {
      loop <- dplyr::filter(loop, s_y >= start_y, e_y <= stop_y, s_x >= start_x, e_x <= stop_x)} else {
        loop <- dplyr::filter(loop, c_y == tad.chr_y, c_x == tad.chr_x, s_y >= start_y, e_y <= stop_y, s_x >= start_x, e_x <= stop_x)}


    p = p +
      ggplot2::geom_rect(data = loop, ggplot2::aes(xmin = s_x, xmax = e_x, ymin = -s_y, ymax = -e_y),
                         fill = NA, color = annotations.color, size = 0.3)
  }

  p
}





















