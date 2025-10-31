#' @title Histogram of annotation distribution within domains
#'
#' @description Plot of the distribution of the genomic annotations within the domains (relative distance)
#' Measures the relative distances of annotation border (start, end) within the domain and plot the distribution as an histogram.
#'
#'
#' @details `domainHist()` measures the relative positions of annotation features (start, end, closest or the center of annotations, more details in `areahist()`) within domains and plot the distribution as an histogram.
#' It is possible that some annotations overlap a boundary, in this cases few option are available to get the relative positions (see `ifoverlap` parameter and the example for better understanding):
#'     * remove those annotations,
#'     * uses the domain in which the feature (`"start"`, `"end"` or `"center"` of the annotation) is located,
#'     * uses the domain in which the annotation has the best overlap.
#'     Therefore, in some cases, the `"start"` of an `annot.gr` can be located before the domain (i.e before the domain with the best overlap). In that case, the distance between the `"start"` of `annot.boundary` and the domain is represented with the real distance (in base pair).
#'
#'
#' @inheritParams boundArea
#' @param annot.boundary Type of feature to analyzed. `"Start"`, `"end"` or `"center"` of each `annot.gr`. It is possible to use `"closest"` to take the closest border (start or end) independently to the strands.
#' @param annot.strand Default is `FALSE` to plot the distribution as histogram. If `TRUE`, distributions are separated according to their strands and are displayed with lines.
#' @param bin.width Size of the bin in percent to count the number of annotations features, default is 5%. if `ifoverlap = "best"`: the real bin distances (measured outside of domains) are equal to `bin.width * 1000`. Therefore default real bin size is 5kb.
#' @param ifoverlap In case of annotation overlap a TAD boundary, few options are available to measure the `annot.boundary` positions:
#'    * `"remove"` to remove all `annot.gr` that overlaps a TAD boundary,
#'    * `"real"` to take the position according to the TAD where the `annot.boundary` is located,
#'    * `"best"` to take the position according to the TAD where the `annot.gr` has the best overlay.
#' @param output Default is `"plot"` to return a `ggplot`. Use `"data"` to return the datas used to produce the plot.
#'
#' @return `ggplot` object.
#'
#' @import GenomeInfoDb
#' @importFrom plyr ddply
#' @importFrom stats quantile
#' @import ggplot2
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' # see vignette("Turorial_TADkit_R_package") or on github (https://github.com/Nico-FR/TADkit)
#'
#' #Create 20 genes
#' annot.gr <- dataframes2grange(
#'   data.frame(chr = 1, start = seq(190.5e3, 209.5e3, 1e3), end = 290e3, strand = "+"),
#'   data.frame(chr = "1", size = 400e3),
#'   strand.col = 4
#' )
#' #Create 2 TADs
#' tad.gr <- dataframes2grange(
#'   data.frame(chr = 1, start = c(100e3, 200e3), end = c(200e3, 300e3)),
#'   data.frame(chr = "1", size = 400e3)
#' )
#'
#' TADplot(tad.gr = tad.gr, annot.gr = annot.gr, start = 150e3, stop = 300e3, chr = 1)
#'
#' #see distribution of gene starts according to the "ifoverlap" parameters:
#' domainHist(domain.gr = tad.gr, annot.gr = annot.gr, ifoverlap = "real")
#' domainHist(domain.gr = tad.gr, annot.gr = annot.gr, ifoverlap = "remove")
#' domainHist(domain.gr = tad.gr, annot.gr = annot.gr, ifoverlap = "best")
#'
domainHist <- function(domain.gr, annot.gr, annot.boundary = "start", ifoverlap = "remove", annot.strand = FALSE, bin.width = 5, output = "plot") {

  if (is.na(mean(seqlengths(domain.gr), na.rm=T))) {
    stop("domain.gr must have seqlengths datas (see dataframes2grange function)")
  }

  # create interTAD annotation and merged them with domain.gr
  interdomain.gr <- GenomicRanges::gaps(domain.gr)[BiocGenerics::strand(GenomicRanges::gaps(domain.gr)) == "*"]
  interdomain.gr$TAD <- FALSE
  domain.gr$TAD <- TRUE
  genome.gr <- c(interdomain.gr, domain.gr)


  annot.gr$nb_overlap_tad <- GenomicRanges::countOverlaps(annot.gr, domain.gr) # number of TAD overlapped
  annot.gr$nb_overlap_gap <- GenomicRanges::countOverlaps(annot.gr, interdomain.gr) # number of interTAD overlapped
  annot.gr$nb_overlap <- GenomicRanges::countOverlaps(annot.gr, genome.gr) # number of TAD or interTAD overlapped (ie boundaries overlapped)

  message(paste0(length(annot.gr[annot.gr$nb_overlap_tad == 0]),"/", length(annot.gr)," annotations are outside domains"))
  message(paste0(length(annot.gr[annot.gr$nb_overlap_tad >= 2]) + length(annot.gr[annot.gr$nb_overlap_tad == 1 & annot.gr$nb_overlap_gap >= 1]), "/", length(annot.gr)," annotations are overlapping with a boundary"))
  message(paste0(length(annot.gr[annot.gr$nb_overlap_tad == 1 & annot.gr$nb_overlap_gap == 0]), "/", length(annot.gr)," annotations are within domains and do not overlap a boundary"))

  annot2.gr <- GenomicRanges::resize(annot.gr, 1, fix = annot.boundary) # keep annot.boundary

  # add (inter)TADhit (TAD line number) in which the annot.boundary is in:
  annot2.gr$genomeHit <- GenomicRanges::findOverlaps(annot2.gr, genome.gr, select = "arbitrary") # arbitrary to randomly select the (inter)TAD when annot.boundary have the same position than a boundary

  annot2.gr$TADhit <- genome.gr$TAD[annot2.gr$genomeHit] # annot.boundary in TAD or interTAD?


  ###################################
  # real: distribution according to annot.boundary positions of each TADs
  ###################################
  if (ifoverlap == "real") {
    data.gr <- annot2.gr[!is.na(annot2.gr$TADhit) & annot2.gr$TADhit == "TRUE"] # remove genes which are not overlapping a TAD
  }
  ###################################

  ###################################
  # remove: distribution according to annot.boundary position without annot.gr that overlap a boundary (2+ TADs or 1+ TADs & 1+ gaps)
  ###################################
  if (ifoverlap == "remove") {
    data.gr <- annot2.gr[annot2.gr$nb_overlap_tad == 1 & annot2.gr$nb_overlap_gap == 0]
  }
  ###################################

  ###################################
  # best: distribution according to annot.boundary position using the best match of annot.gr
  ###################################
  if (ifoverlap == "best") {
    annot.overlap <- annot.gr[annot.gr$nb_overlap >= 2] # annot.gr (with >=2 overlap)

    # overlap between annot.gr and genome.gr
    over <- as.data.frame(GenomicRanges::findOverlaps(annot.overlap, genome.gr))

    # size of each overlap (annot.overlap VS genome.gr)
    over$over_size <- IRanges::width(GenomicRanges::pintersect(
      annot.overlap[over$queryHits],
      genome.gr[over$subjectHits]
    ))

    # best match
    best <- plyr::ddply(over, ~queryHits, function(d) d[which.max(d$over_size), ])[, 1:2]
    annot.overlap$genomeHit <- best$subjectHits
    annot.overlap$TADhit <- genome.gr$TAD[best$subjectHits] # best match is TAD or interTAD?

    # merged annot2.gr with only one overlap with the best overlap for the other one (+resize according to annot.boundary)
    data.gr <- c(
      annot2.gr[annot2.gr$nb_overlap_tad == 1 & annot2.gr$nb_overlap_gap == 0],
      GenomicRanges::resize(annot.overlap[annot.overlap$TADhit], 1, fix = annot.boundary)
    )
  }
  ###################################



  data2.gr <- data.gr # keep annotations in TADs


  data2.gr$startHit <- GenomicRanges::start(genome.gr[data2.gr$genomeHit]) # (inter)TAD start
  data2.gr$widthHit <- BiocGenerics::width(genome.gr[data2.gr$genomeHit]) # (inter)TAD width

  # relative position (%) of annot.boundary according to (inter)TAD start and (inter)TAD size
  data2.gr$relative_position <- (start(data2.gr) - data2.gr$startHit) / data2.gr$widthHit * 100

  # real position of annot.boundary outside of (inter)TAD (otherwise = NA)
  data2.gr$real_position <- ifelse(start(data2.gr) < data2.gr$startHit,
    start(data2.gr) - data2.gr$startHit,
    ifelse(start(data2.gr) > (data2.gr$startHit + data2.gr$widthHit), start(data2.gr) - data2.gr$startHit - data2.gr$widthHit, NA)
  )

  # mixed position: real:[-Inf, 0[ ; relative:[0,100] ; real:]100, +Inf]
  data2.gr$mixed_position <- ifelse(data2.gr$relative_position < 0, data2.gr$real_position,
    ifelse(data2.gr$relative_position > 100, data2.gr$real_position + 100 * 1e3, data2.gr$relative_position * 1e3)
  )

  ###############################################
  #OUTPUTS
  ###############################################
  if (output == "data") {

    data2.gr$nameHit = names(genome.gr)[data2.gr$genomeHit]
    data2.gr$genomeHit = NULL
    data2.gr$TADhit = NULL
    return(BiocGenerics::sort(data2.gr, ignore.strand = TRUE))

  }
  if (output == "plot") {

  ##############
  # histo
  ##############
  names(data2.gr) <- NULL # some names could be duplicated, so remove the names to create a data.frame

  if (isFALSE(annot.strand)) {
    p <- ggplot2::ggplot(as.data.frame(data2.gr), ggplot2::aes(x = mixed_position)) +
      ggplot2::geom_histogram(binwidth = bin.width * 1e3, boundary = 0, size = 0.75, color = "black", fill = "white")
  }

  if (isTRUE(annot.strand)) {
    p <- ggplot2::ggplot(as.data.frame(data2.gr), ggplot2::aes(x = mixed_position, color = strand)) +
      ggplot2::geom_freqpoly(binwidth = bin.width * 1e3, boundary = 0, size = 0.75, alpha = 0.75) +
      ggplot2::scale_color_manual(values = c("#33A02C", "#1F78B4")) +
      ggplot2::geom_point(stat = "bin", ggplot2::aes(y = after_stat(count)), binwidth = bin.width * 1e3, boundary = 0, size = 1, alpha = 0.75)
  }


  # number of 20kb intervals upstream
  n_end_limit <- ceiling((quantile(data2.gr$mixed_position, na.rm = T, probs = 0.99) - 100e3) / 20e3)

  # number of 20kb intervals downstream
  n_start_limit <- -1 * floor(quantile(data2.gr$mixed_position, na.rm = T, probs = 0.01) / 20e3)


  return(
    p + ggplot2::geom_vline(ggplot2::aes(xintercept = 0), color = "red", size = 0.75, alpha = 0.75) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = 100 * 1e3), color = "red", size = 0.75, alpha = 0.75) +
      ggplot2::ggtitle(paste0(annot.boundary, "::", ifoverlap)) + ggplot2::xlab("") +
      ggplot2::scale_x_continuous(
        breaks = c(seq(
          ifelse(n_start_limit > 0, -n_start_limit * 20e3, 0),
          ifelse(n_end_limit > 0, 100e3 + n_end_limit * 20e3, 100e3), 20e3
        )),
        limits = c(
          ifelse(n_start_limit > 0, -n_start_limit * 20e3, 0),
          ifelse(n_end_limit > 0, 100e3 + n_end_limit * 20e3, 100e3)
        ),
        labels = c(
          if (n_start_limit > 0) {
            paste0((1:n_start_limit * 20) * -1, "kb")
          },
          paste0(seq(0, 100, 20), "%"),
          if (n_end_limit > 0) {
            paste0("+", (1:n_end_limit * 20), "kb")
          }
        )
      )
  )
}}






