#' plot view point interaction
#'
#' @description This function calculates the average number of interactions (i.e. expected number of interaction) as a function of the distance between bins.
#'
#' @inheritParams mMATplot
#' @param matrix.lst List of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param start,stop Region of interest in base pair. Default is NULL to use the entire chromosome (i.e. entire matrix).
#' @param vp.start Start of the view point in base pair.
#' @param vp.stop Stop/end of the view point in base pair.
#' @param self_interaction logical. Whether or not to add interactions within the view-point.
#' @param norm Normalized the interaction count by the average interaction of the view point (mean interaction along the chromosome). Default is FALSE.
#' @param log2 logical. Use the log2 of the matrix values. Default is `FALSE`. Note that if TRUE, interaction with 0 count are removed from the analysis.
#' @return ggplot
#'
#' @importFrom scales unit_format
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' # get domains from boundaries:
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' domains.gr = boundary2domain(boundaries.gr)
#'
#' MATplot(mat_HCT116_chr19_50kb,
#'     start = 7e6, stop = 15e6,
#'     bin.width = 50e3, log2 = TRUE,
#'     tad.upper.tri = domains.gr["chr19_11597500",])
#'
#' #plot interactions of the domain above, between 10Mb and 15Mb
#' viewPointInteract(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb), bin.width = 50e3,
#'   vp.start = 11597500, vp.stop = 12672500,
#'   start = 7e6, stop = 15e6,
#'   self_interaction = FALSE) #FALSE to not plot intra-domain interactions
#'
#'
#'
#'
viewPointInteract <- function(matrix.lst, bin.width, vp.start, vp.stop, start = NULL, stop = NULL, self_interaction = FALSE, norm = FALSE, log2 = FALSE,
                       colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")) {

  #sanity check
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  #matrix coordinate of limits
  if (is.null(start) & is.null(stop)) {
    stop = nrow(matrix.lst[[1]])} else {
      l.start = start %/% bin.width + 1
      l.stop = stop %/% bin.width
    }

  #matrix coordinate of view point
  if ((vp.start %/% bin.width == vp.start / bin.width) & (vp.stop %/% bin.width == vp.stop / bin.width)) {
    vp.start = vp.start / bin.width
    vp.stop = vp.stop / bin.width} else {
      vp.start = round(vp.start / bin.width) + 1 # +1 because matrix coordinates start to 1. So first nucleotide (i.e. vp.start = 1bp) is bin nb 1.
      vp.stop = ifelse((round(vp.stop / bin.width) + 1) <= vp.start, vp.start, round(vp.stop / bin.width))
      message(paste0("vp.start/vp.stop are not multiples of bin.width, round to ", (vp.start - 1) * bin.width,
                     " and ", vp.stop * bin.width, " (i.e. ", vp.stop + 1 - vp.start, " bins)."))
    }

  #upper matrix to symetrical matrix (need to be improve by skeeping symetrical matrix)
  matrix2.lst = lapply(matrix.lst, function(m) {
    mat = as.matrix(m)
    mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
    #diag(mat) <- NA
    return(mat[vp.start:vp.stop,])
    })

  #norm = TRUE
  if (isTRUE(norm)) {
    #observed / expected (expected = mean interaction of all view point bins)
    matrix2.lst = mapply("/", matrix2.lst, mapply(base::mean, matrix2.lst, SIMPLIFY = FALSE), SIMPLIFY = FALSE)
  }

  #log2 = TRUE: remove 0 counts
  if (isTRUE(log2)) {
    matrix2.lst = lapply(matrix2.lst, function(m) {
      m[m == 0] <- NA
      return(m)}) #remove 0 count
    matrix2.lst = lapply(matrix2.lst, function(m) {log2(m)})
  }

  #colwise mean
  mean_interact.lst = if (vp.start == vp.stop) {
    lapply(matrix2.lst, function(m) {m[l.start:l.stop]})
    } else {
      lapply(matrix2.lst, function(m) {colMeans(m[,l.start:l.stop], na.rm = TRUE)})
      }

  #list to data frame
  df = do.call(data.frame, mean_interact.lst)

  #self_interaction = TRUE
  if (isFALSE(self_interaction)) {
    #remove view point self interaction
    df[(vp.start - l.start + 1):(vp.stop - l.start + 1), ] <- NA
  }

  df[nrow(df) + 1,] = df[nrow(df),] #copy last bin to visualized the last bin in the graph (geom_step)
  df$distance = (l.start:(l.stop + 1) - 1) * bin.width #bin start



  if (is.null(names(matrix.lst))) {names(df) = c(1:length(matrix.lst), "distance")}

  p = ggplot2::ggplot(data = tidyr::gather(df, "matrix",  "vp_interaction_avg", 1:length(matrix.lst)),
                      aes(y = vp_interaction_avg, x = distance, color = matrix))+
    ggplot2::geom_step()+scale_color_manual(values = colors.lst)+
    scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
    geom_vline(xintercept = c((vp.start - 1) * bin.width, vp.stop * bin.width), linetype = "dotted")+
    ylab(ifelse(isTRUE(log2), "log2(vp_interaction_avg)", "vp_interaction_avg"))

  return(p)
}

















