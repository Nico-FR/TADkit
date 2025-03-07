#' plot view point interaction
#'
#' @description This function (also known as virtual4C) calculates the average number of interactions of a region of interest (e.g. a TAD) with respect to the rest of the matrix.
#'
#' @inheritParams mMATplot
#' @param matrix.lst List of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param output "GRanges" (default) to return a bedgraph with the average interaction of the view point with the area analyzed and "plot" to return a graph.
#' @param vp.start Start of the view point in base pair.
#' @param vp.stop Stop/end of the view point in base pair.
#' @param self_interaction logical. Whether or not to add interactions within the view-point.
#' @param start,stop Area in bp of the chromosome to be analyzed. Default is NULL to use the entire chromosome (i.e. entire matrix).
#' @param Qnorm Quantile normalization of the average interaction. Default is FALSE.
#' This parameter can be useful for the graph, as it allows a better comparison of distributions between matrices. Note that normalization is best performed on the entire chromosome than anly on an area.
#' @param log2 logical. Use the log2 of the matrix values. Default is `FALSE`. Note that if TRUE, interaction with 0 count are removed from the analysis.
#' @param colors.lst Set of 8 colors used for plot.
#' @param seqname chromosome names as character, default = "1".
#' @return `GRanges` bedgraph with average interaction.
#'
#' @importFrom scales unit_format
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' # get domains from boundaries:
#' tad_boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' tad_domains.gr = boundary2domain(tad_boundaries.gr)
#'
#' MATplot(mat_HCT116_chr19_50kb,
#'     start = 7e6, stop = 15e6,
#'     bin.width = 50e3, log2 = TRUE,
#'     tad.upper.tri = tad_domains.gr["chr19_11597500",])
#'
#' #plot interactions of the domain above, between 7Mb and 15Mb
#' viewPointInteract(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb), bin.width = 50e3,
#'   vp.start = 11597500, vp.stop = 12672500,
#'   start = 7e6, stop = 15e6,
#'   self_interaction = FALSE, #FALSE to not plot intra-TAD interactions
#'   output = "plot")
#'
#'
viewPointInteract <- function(matrix.lst, bin.width, vp.start, vp.stop, output = "GRanges", start = NULL, stop = NULL, self_interaction = FALSE, Qnorm = FALSE, log2 = FALSE,
                              seqname = "1", colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")) {

  #matrix.lst = mat.lst
  #bin.width = 4e3
  #vp.start = 2380000; vp.stop = 2730000
  #output = "plot"; self_interaction = F; Qnorm = F;log2 = TRUE
  #start = 2.06e6; stop = 3.06e6

  mean_interact <- . <- NULL
  #sanity check
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  if (is.null(names(matrix.lst))) {names(matrix.lst) = 1:length(matrix.lst)}

  #matrix coordinate of limits
  if (is.null(start)) {
    l.start = 1}
  if (is.null(stop)) {
    l.stop = nrow(matrix.lst[[1]])}

  if (!is.null(start)) {
    l.start = start %/% bin.width + 1}
  if (!is.null(stop)) {
    l.stop = stop %/% bin.width}

  #matrix coordinate of view point
  if ((vp.start %/% bin.width == vp.start / bin.width) & (vp.stop %/% bin.width == vp.stop / bin.width)) {
    vp.start = vp.start / bin.width + 1 #+1 because bin 1 start to 0 and stop to bin.width
    vp.stop = vp.stop / bin.width
  } else {
    vp.start = round(vp.start / bin.width) + 1 # +1 because matrix coordinates start to 1. So first nucleotide (i.e. vp.start = 1bp) is bin nb 1.
    vp.stop = ifelse((round(vp.stop / bin.width)) <= vp.start, vp.start, round(vp.stop / bin.width))
    message(paste0("vp.start/vp.stop are not multiples of bin.width, round to ", (vp.start - 1) * bin.width,
                   " and ", vp.stop * bin.width, " (i.e. ", vp.stop + 1 - vp.start, " bins)."))
  }

  #crop matrix to vp
  matrix2.lst = lapply(matrix.lst, function(MAT) {

    mat = as.matrix(MAT[l.start:l.stop, l.start:l.stop])

    if (!isSymmetric(mat)) {
      mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    }

    return(mat[(vp.start - l.start + 1):(vp.stop - l.start + 1),])})

  #log2 = TRUE: remove 0 counts
  if (isTRUE(log2)) {
    matrix2.lst = lapply(matrix2.lst, function(MAT) {
      MAT[MAT == 0] <- NA
      return(log2(MAT))})}

  #colwise mean
  mean_interact.mat = if (vp.start == vp.stop) {
    sapply(matrix2.lst, function(m) {m})
  } else {
    sapply(matrix2.lst, function(m) {colMeans(m, na.rm = TRUE)})
  }

  #option
  #quantile normalization
  if (Qnorm == TRUE) {
    preprocessCore::normalize.quantiles(mean_interact.mat, copy=FALSE, keep.names=TRUE)
  }
  #self_interaction
  if (isFALSE(self_interaction)) {
    #remove view point self interaction
    mean_interact.mat[(vp.start - l.start + 1):(vp.stop - l.start + 1), ] <- NA
  }

  mean_interact.df = data.frame(seqname = seqname,
                                start = (l.start:l.stop - 1) * bin.width,
                                stop = l.start:l.stop * bin.width) %>%
    cbind.data.frame(as.data.frame(mean_interact.mat))

  df.lst = lapply(1:ncol(mean_interact.mat), function(i) {
    df = data.frame(seqname = seqname,
                    start = (l.start:l.stop - 1) * bin.width,
                    stop = l.start:l.stop * bin.width,
                    mean_interact = mean_interact.mat[,i])
  })
  names(df.lst) = colnames(mean_interact.mat)

  colname = if(isTRUE(log2) & isFALSE(Qnorm)){"log2_mean_interact"
  } else if(isTRUE(log2) & isTRUE(Qnorm)){"log2_Qnorm_mean_interact"
  } else if(isFALSE(log2) & isTRUE(Qnorm)){"Qnorm_mean_interact"
  } else if(isFALSE(log2) & isFALSE(Qnorm)){"mean_interact"}

  if (output == "GRanges") {
    gr.lst = lapply(df.lst, function(DF) {
      DF %>% `colnames<-`(c(names(df.lst[[1]][1:3]), colname)) %>%
        GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)})
    return(gr.lst)
  }

  if (output == "plot") {
    lapply(1:length(df.lst), function(INT) {
      df.lst[[INT]] %>%
        dplyr::add_row(seqname = seqname,
                       start = df.lst[[INT]][nrow(df.lst[[INT]]), 3],
                       stop = df.lst[[INT]][nrow(df.lst[[INT]]), 3],
                       mean_interact = df.lst[[INT]][nrow(df.lst[[INT]]), 4]) %>% #copy last row to visualized the last bin in the graph (geom_step)
        dplyr::mutate(matrix = names(df.lst[INT])) %>%
        dplyr::rename("distance" = "start")
    }) %>%
      do.call(rbind, .) %>%
      ggplot2::ggplot(., aes(y = mean_interact, x = distance, color = matrix))+
      ggplot2::geom_step()+scale_color_manual(values = colors.lst)+
      scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
      geom_vline(xintercept = c((vp.start - 1) * bin.width, vp.stop * bin.width), linetype = "dotted")+
      ylab(colname) -> p
    return(p)
  }
}



