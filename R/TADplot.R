#' @title Plot domains for one individual
#'
#' @description `TADplot()` use `Gviz` package to plot at least 2 tracks:
#' * All the TADs from a specified chromosome,
#' * A zoom in the defined area.
#'
#' Other tracks can be added:
#' * bed files with annotations (genes...),
#' * bigwig files (read depth sequencing...),
#' * bedgraph with bin values (insulation score...).
#'
#' @details This create a plot with at least 2 tracks from TADs annotation (tad.gr).
#'
#' On the first track the alternating TADs are represented in grey and black while the interTAD regions are represented in white.
#' The second track is a zoom in the area of interest. This area is extend to the border of the TADs on both sides.
#'
#' Another track from a bigwig file can be added like read depth sequencing (RNAseq...) as an histogram. However bigwig files can not be read on Windows.
#' The bin size of the histogram is 1Kb by default, the read depth is smoothed using the median value of each bin.
#' If the chromosome name is different (between the bigwig and the `chr` parameter), it can be fix using the `bigwig.chr` parameter (e.g `chr="1"` while `bigwig.chr="chr1"`).
#'
#' Another track with any annotations can be added.
#' This track can be group using factors in a specified column of the `annot.gr` files (metadata) otherwise the names of each annotation is used.
#'
#' @param tad.gr `GRanges` object with domains. The chromosomes lengths must be stored as metadata (see `dataframes2grange()`).
#' @param chr Chromosome name to plot.
#' @param start,stop Region of interest in base pair.
#' @param tad.id Logical. Default is `FALSE` to label the domain with their sizes `width(tad.gr)` instead of their names `names(tad.gr)`.
#' @param bigwig.path Path for the bigwig file plotted as histogram. Default = `NULL` (i.e the track is not plotted). Note that bigwig file can not be read on Windows.
#' @param annot.gr `GRanges` object with genomic annotations. Default = `NULL` (i.e the track is not plotted).
#' @param bigwig.binwidth Bin sizes for the histogram of the bigwig track. Default = `1e3`.
#' @param bigwig.xaxis Function used to transform the x-axis of the bigwig values among each bigwig.binwidth. Defaults = `"median"`.
#' Alternatively, other predefined functions can be supplied as character (`"mean"`, `"median"`, `"sum"`, `"min"`, `"max"` or `"extreme"`).
#' @param bigwig.chr Chromosome name used for the bigwig file. Default = `NULL` to used the same name as chr.
#' @param bigwig.yaxis Function used to transforming the y-axis of bigwig values. Default = `NULL`. Use `"log2"` to use the function `log2(x + 1)` to transform the y-axis or provide any other function.
#' @param annot.col Column number of the metadata from `annot.gr` file used to group the annotation tracks. Default = `NULL`.
#' @param bedgraph `data.frame`, `GRanges` or path for the bedgraph file (data frame without header and 4 columns tab separated) plotted as line. Default = NULL (i.e no track is plotted).
#'
#' @return Plot with domains and other tracks as a list of GenomeGraph tracks (see `Gviz::plotTracks` for details).
#'
#' @import GenomicRanges
#' @importFrom S4Vectors aggregate
#' @importFrom IRanges IRanges
#' @import GenomeInfoDb
#' @import Gviz
#' @import rtracklayer
#'
#' @export
#'
#' @examples # Lets create a GRange object with 3 TADs and 1 metadata column:
#' tad.gr <- dataframes2grange(
#'   data.frame(
#'     chr = rep(1, 3),
#'     start = c(2e6, 3e6, 4e6),
#'     end = c(3e6, 4e6, 5e6),
#'     type = c("a", "b", "a")
#'   ),
#'   data.frame(chr = "1", size = 6e6),
#'   metadata.mcols = 4
#' )
#' tad.gr
#'
#' # plot those TADs
#' # use those TADs as annotation tracks (with the metadata column to group annotations):
#' TADplot(
#'   tad.gr = tad.gr,
#'   chr = 1, start = 2e6, stop = 5e6,
#'   annot.gr = tad.gr, annot.col = 1
#' )
TADplot <- function(tad.gr, chr, start, stop, tad.id = FALSE,
                    bigwig.path = NULL, bigwig.binwidth = 1e3, bigwig.xaxis = "mean", bigwig.chr = NULL, bigwig.yaxis = NULL,
                    annot.gr = NULL, annot.col = NULL, bedgraph = NULL) {

  options(ucscChromosomeNames=FALSE)

  ##############################
  # Ideotrack
  ##############################
  data <- tad.gr[GenomeInfoDb::seqnames(tad.gr) == chr]
  data$gieStain <- ifelse(1:length(data) %% 2, "gpos100", "gpos25") # alternate colors for TAD
  gaps <- GenomicRanges::gaps(data)
  gaps <- gaps[GenomeInfoDb::seqnames(gaps) == chr & BiocGenerics::strand(gaps) == "*"]
  if (length(gaps) != 0) {gaps$gieStain <- "gpos1"} # color for interTAD


  data2 <- sort(c(data, gaps))
  d1 <- data.frame(
    chrom = paste0("chr", gsub('chr','', chr)),
    chromStart = GenomicRanges::start(data2),
    chromEnd = BiocGenerics::end(data2),
    name = paste0("n", 1:length(data2[,1])),
    gieStain = data2$gieStain
  )

  ideoTrack <- Gviz::IdeogramTrack(genome = "custom", chromosome = paste0("chr", gsub('chr','', chr)), fontsize = 10, bands = d1)

  ##############################
  # tadTrack
  ##############################
  data <- tad.gr[GenomeInfoDb::seqnames(tad.gr) == chr &
                   BiocGenerics::end(tad.gr) >= start &
                   GenomicRanges::start(tad.gr) <= stop]

  if (length(data)==0) {tadTrack = Gviz::GenomeAxisTrack(
    add53 = T,
    add35 = T,
    littleTicks = T,
    lwd = 1,
    col = "black"
  )
  extended_start <- start
  extended_stop <- stop
  } else {
    extended_start <- min(start(data), start)
    extended_stop <- max(end(data), stop)

    tadTrack <- Gviz::GenomeAxisTrack(
      range = IRanges::IRanges(
        start = GenomicRanges::start(data) + 1,
        end = BiocGenerics::end(data) - 1,
        names = if (!isTRUE(tad.id)) {paste0(round(BiocGenerics::width(data) / 1e3), "Kb")} else {names(data)}
      ),
      fill.range = c("cornflowerblue", "coral1"), add53 = T,
      add35 = T, littleTicks = T,
      lwd = 1, col = "black", showId = T, col.id = "black"
    )
  }

  ##############################
  # bigwigTrack
  ##############################
  if (.Platform$OS.type == "windows") {bigwig.path = NULL}
  if (is.null(bigwig.path)) {
    bigwigTrack <- NULL
  } else {

    # chr name to use
    if (is.null(bigwig.chr)) {
      bigwig.chr <- chr
    }

    if (is.null(bigwig.yaxis)) {
      transformation.method = NULL} else {
        if (deparse(bigwig.yaxis)[1] == deparse("log2")[1]) {
          transformation.method = function(x) { log2(x + 1) }
        } else {transformation.method = bigwig.yaxis}
      }

    # creation of the import region as GRange
    bin.gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = extended_start, end = extended_stop),
      strand = "*"
    )
    GenomeInfoDb::seqlengths(bin.gr) <- GenomeInfoDb::seqlengths(tad.gr)[as.character(chr)]

    GenomeInfoDb::seqlevels(bin.gr) <- as.character(bigwig.chr)

    # import bigwig files as Rle
    cov.rle <- rtracklayer::import.bw(bigwig.path,
                                      which = bin.gr, as = "RleList"
    )

    # Sanity check of chr names of bigwig file
    if (length(cov.rle[names(cov.rle) == as.character(bigwig.chr)]) == 0) {
      stop(paste0("There is no chromosome ", bigwig.chr, " in the bigwig file."))
    }

    bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(bin.gr)[as.character(bigwig.chr)],
                                      tilewidth = bigwig.binwidth, cut.last.tile.in.chrom = TRUE
    )

    bins2 <- bins[start(bins) >= extended_start & end(bins) <= extended_stop]

    # coverage for each bins according to aggregation.method
    temp <- S4Vectors::aggregate(cov.rle[[as.character(bigwig.chr)]], bins2, FUN = bigwig.xaxis)
    bins2$median_cov <- temp # create metadata with the median values

    GenomeInfoDb::seqlevels(bins2) = paste0("chr", gsub('chr','', GenomeInfoDb::seqlevels(bins2)))

    bigwigTrack <- Gviz::DataTrack(bins2,
                                   chr =  paste0("chr", gsub('chr','', chr)),
                                   type = "hist",
                                   fill = "grey50",
                                   name = "bw",
                                   col.histogram = "grey25",
                                   transformation = transformation.method
    )
  }

  #####################################
  # annotTracks
  #####################################
  if (is.null(annot.gr)) {
    annotTrack <- NULL
  } else {
    if (is.null(annot.col)) {
      data1 <- annot.gr[GenomeInfoDb::seqnames(annot.gr) == chr &
                          BiocGenerics::end(annot.gr) >= extended_start &
                          GenomicRanges::start(annot.gr) <= extended_stop]

      if (length(data1) == 0) {
        annotTrack <- NULL
      } else {
        annotTrack <- Gviz::AnnotationTrack(
          start = GenomicRanges::start(data1),
          width = BiocGenerics::width(data1),
          chromosome = paste0("chr", gsub('chr','', chr)),
          strand = BiocGenerics::strand(data1),
          id = names(data1), groupAnnotation = "id",
          genome = "custom", name = "Annot", col = "deepskyblue4"
        )
      }
    } else {
      annotTrack <- NULL

      # Sanity check of annot.col
      if (length(GenomicRanges::mcols(annot.gr)) < annot.col) {
        stop(paste0("Wrong annot.col number. There is only ", length(GenomicRanges::mcols(annot.gr)), " column(s) with metadata."))
      }

      data1 <- annot.gr[GenomeInfoDb::seqnames(annot.gr) == chr &
                          BiocGenerics::end(annot.gr) >= extended_start &
                          GenomicRanges::start(annot.gr) <= extended_stop]

      if (length(data1) == 0) {
        annotTrack <- NULL
      } else {
        annotTrack <- Gviz::AnnotationTrack(
          start = GenomicRanges::start(data1),
          width = BiocGenerics::width(data1),
          chromosome = paste0("chr", gsub('chr','', chr)),
          strand = BiocGenerics::strand(data1),
          group = GenomicRanges::mcols(data1)[,annot.col], groupAnnotation = "group",
          genome = "custom", name = "Annot", col = "deepskyblue4",
          col.line = "darkblue", cex.feature = 0.5, cex.group = 0.7,
          just.group = "below"
        )
      }
    }
  }

  #####################################
  # bedgraphTrack
  #####################################
  if (is.null(bedgraph)) {
    bedgraphTrack <- NULL
  } else {
    chromsize = data.frame(chr = GenomeInfoDb::seqlevels(tad.gr),
                           size = GenomeInfoDb::seqlengths(tad.gr))
    #read bedgraph datas
    ##if dataframe
    if (is.data.frame(bedgraph)) {
      data1 = TADkit::dataframes2grange(bedgraph, chromsize, metadata.mcols = 4)
    }
    ##if path
    if (is.character(bedgraph)) {
      data0 = utils::read.table(bedgraph, header = FALSE, sep = "\t")[,1:4]
      data1 = TADkit::dataframes2grange(data0, chromsize, metadata.mcols = 4)
    }
    ##if GRanges
    if (inherits(bedgraph, "GRanges")) {
      data1 = bedgraph[,1] %>% GenomeInfoDb::keepSeqlevels(chr, pruning.mode = "coarse")
    }

    GenomeInfoDb::seqlevels(data1) = paste0("chr", gsub('chr','', GenomeInfoDb::seqlevels(data1)))

    bedgraphTrack <- DataTrack(range = data1, type = "b", lwd = 1, name = "bg", chromosome =  paste0("chr", gsub('chr','', chr)))
  }

  #####################################
  # grid
  #####################################
  data = GenomicRanges::restrict(tad.gr[GenomeInfoDb::seqnames(tad.gr) == chr],
                                 start = extended_start, end = extended_stop)

  ht <- Gviz::HighlightTrack(trackList = c(tadTrack, bedgraphTrack, bigwigTrack, annotTrack),
                             start = unique(sort(c(GenomicRanges::start(data),extended_start, BiocGenerics::end(data),extended_stop))),
                             end = unique(sort(c(GenomicRanges::start(data),extended_start, BiocGenerics::end(data),extended_stop))),
                             chr = paste0("chr", gsub('chr','', chr)), col = "grey70", fill = "white")

  #####################################
  # Plot
  #####################################
  Gviz::plotTracks(c(ideoTrack, ht),
                   from = extended_start,
                   to = extended_stop,
                   background.title = "grey30", grid = TRUE, v = 0
  )
}
