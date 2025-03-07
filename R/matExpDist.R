#' expected interactions according to genetic distances
#'
#' @description This function calculates the average number of interactions (i.e. expected number of interaction) as a function of the distance between bins.
#'
#' @param matrix.lst List of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Bin width of the matrices in base pair.
#' @param colors.lst Set of 8 colors used for plot.
#' @return ggplot
#'
#' @importFrom scales unit_format
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom Matrix summary
#' @importFrom tidyr complete
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' matExpDist(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb), bin.width = 50e3)
#'
#'
matExpDist <- function(matrix.lst, bin.width,
                       colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")) {

  i <- j <- expected <- distance <- matrix <- . <- x <- NULL
  #sanity check
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  matrix.lst = lapply(matrix.lst, function(MAT) {
    if(inherits(MAT, "matrix")) {
      methods::as(MAT, "CsparseMatrix")} else {MAT}
  })

  output = lapply(1:length(matrix.lst), function(INT) {
    Matrix::summary(matrix.lst[[INT]]) %>%
      dplyr::mutate(i = factor(i, levels = 1:ncol(matrix.lst[[INT]]))) %>% #add missing index as levels
      dplyr::mutate(j = factor(j, levels = 1:ncol(matrix.lst[[INT]]))) %>% #add missing index as levels
      tidyr::complete(., i, j, fill = list(x = NA), explicit = FALSE) %>% #add NA values to missing bins
      dplyr::filter(as.numeric(i) <= as.numeric(j)) %>% #filter lower matrix
      dplyr::mutate(distance = abs(as.numeric(i) - as.numeric(j)) * bin.width) %>% #add distances between bins
      dplyr::group_by(distance) %>% dplyr::summarise(expected = mean(x, na.rm = TRUE)) %>% # mean according to distances
      dplyr::filter(!is.na(expected)) %>%
      dplyr::mutate(matrix = ifelse(is.null(names(matrix.lst[INT])), paste0("matrix", INT), names(matrix.lst[INT])))
  })

  df = do.call(rbind, output)

  p = ggplot2::ggplot(data = df, ggplot2::aes(y = expected, x = distance, color = matrix))+
    ggplot2::geom_line()+ggplot2::scale_color_manual(values = colors.lst)+
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), trans = "log10")+
    ggplot2::scale_y_continuous(trans = "log10")

  return(p)
}






































