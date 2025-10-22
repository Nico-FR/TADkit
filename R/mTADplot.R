#' Plot domains for multiple individuals
#'
#' @description `mTADplot()` is similar to `TADplot()` function but it allows to plot and compare several individuals.
#' The simplest graph allows to view the domains of one or more individuals.
#'
#' Other tracks can be added:
#' * bed files with annotations (genes...),
#' * bigwig files (read depth sequencing...),
#' * bedgraph with bin values (insulation score...).
#'
#' @details `mTADplot()` create a plot with at least the track with TADs annotations (from 1 or more files).
#' TAD annotation(s) must be in a list and element in the list must have name(s).
#' TAD annotation(s) must have chromosome size (see `dataframes2grange()`).
#'
#' Another track from a list of bigwig file(s) can be added, like read depth sequencing (RNAseq...) as an histogram. However bigwig files can not be read on Windows.
#' The bin size of the histogram is 1Kb by default, the read depth is smoothed using the median value of each bin.
#' If the chromosome name is different (between the bigwig file(s) and the `chr` parameter), it can be fix using the `bigwig.ch`r parameter (e.g `chr="1"` while `bigwig.chr="chr1"`).
#'
#' Another track from a list of `GRanges` object(s) with any annotations can be added.
#' This track(s) can be group using factors in a specified column of each `GRanges` (metadata), otherwise the names of each annotation is used.
#'
#' Another track from a list of bedgraph file(s) (file of 4 columns, chr, stat, stop and value) can be added. The list must contain the path (not the file) of each bedgraph and each path must have a name.
#' It is also possible to use a list containing another list of bedgraph path. In this case, the first level of list will be represented in a different track.
#'
#' @inheritParams TADplot
#' @param tad.lst List of GRange object with domains. Those files must have chromosomes lengths (see dataframes2grange function).
#' @param bigwigPath.lst List of path for the bigwig file(s) plotted as histogram. Default = NULL (ie no track is plotted). Note that bigwig files can not be read on Windows..
#' @param bigwig.binwidth Bin sizes for the histogram of the bigwig track. Default = 1e3.
#' @param annot.lst List of GRange file(s) with genomic annotations. Default = NULL (ie no track is plotted).
#' @param annot.col Column number of the metadata from annot.gr file(s) used to group the annotation tracks. Default = NULL.
#' @param bedgraph.lst `data.frame`, `GRanges` or list of path for the bedgraph file(s) plotted as line. Default = NULL (i.e no track is plotted).
#' it is possible to create several bedgraph tracks (each containing 1 or more lines) by using a list containing few others lists (see vignette). All list must have names.
#' @param bedgraph.name Name of the bedgraph track when there is only one track (default = "bedgraph"). Otherwise it takes the names of each list.
#' @param bedgraph_outliers Ratio to remove outliers of all bedgraph files. Default is 0 (ie no filter). To remove the first and last percentiles use 0.01.
#' @param colors.lst Set of 8 colors used for each files within a list.
#' @param title Title of the plot. Default is null (ie no title).
#'
#' @return Plot with domains and other tracks as a list of GenomeGraph tracks (see `Gviz::plotTracks` for details).
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges ranges
#' @import GenomeInfoDb
#' @importFrom magrittr %>%
#' @importFrom dplyr select full_join filter mutate
#' @importFrom utils read.table
#' @import Gviz
#' @import rtracklayer
#'
#' @export
#'
#' @examples
#' library(IRanges)
#'
#' # get domains from boundaries:
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' domains.gr = boundary2domain(boundaries.gr)
#'
#' #create list of TADs
#' tad.lst = tad.lst = list(ind1 = domains.gr, ind2 = shift(domains.gr, 50e3))
#'
#' #insulation score (data frame) list
#' IS.gr = dataframes2grange(IS_HCT116_chr19_5kb.bedgraph, human_chromsize, metadata.mcols = 4)
#' IS.lst = list(ind1 = IS.gr, ind2 = shift(IS.gr, 50e3))
#'
#' mTADplot(tad.lst = tad.lst, chr = "chr19",
#'   start = 10e6, stop = 12e6,
#'   bedgraph.lst = IS.lst, bedgraph.name = "IS")
#'
mTADplot <- function(tad.lst, chr, start, stop, tad.id = FALSE,
                     bigwigPath.lst = NULL, bigwig.binwidth = 1e3, bigwig.xaxis = "mean", bigwig.chr = NULL, bigwig.yaxis = NULL,
                     annot.lst = NULL, annot.col = NULL,
                     bedgraph.lst = NULL, bedgraph.name = "bedgraph", bedgraph_outliers = 0,
                     colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"),
                     title = NULL) {

  options(ucscChromosomeNames=FALSE)
  V1 <- V2 <- V3 <- V4 <- . <- NULL

  #sanity check
  if (!is.list(tad.lst)) {
    stop("tad.lst must be a list")
  }

  if (is.na(GenomeInfoDb::seqlengths(tad.lst[[1]])[[as.character(chr)]])) {
    warning("There is no seqlenths in tad.lst, see dataframes2grange function")
    return(GenomeInfoDb::seqlengths(tad.lst[[1]]))
  }

  ##############################
  # Ideotrack
  ##############################
  d1 = data.frame(chrom = paste0("chr", gsub('chr','', names(GenomeInfoDb::seqlengths(tad.lst[[1]])))),
                  chromStart = 0,
                  chromEnd = GenomeInfoDb::seqlengths(tad.lst[[1]]),
                  name = "", gieStain = "gneg")

  ideoTrack <-
    Gviz::IdeogramTrack(genome = "custom", chromosome = paste0("chr", gsub('chr','', chr)),
                        fontsize = 10, bands = d1)
  areaTrack <- Gviz::GenomeAxisTrack(
    add53 = T,
    add35 = T,
    littleTicks = T,
    lwd = 1,
    col = "black"
  )

  ##############################
  # tadTracks
  ##############################
  tadTracks <- NULL
  for (i in 1:length(tad.lst)) {
    data <- tad.lst[[i]][GenomeInfoDb::seqnames(tad.lst[[i]]) == chr]

    Track <- Gviz::AnnotationTrack(
      start = GenomicRanges::start(data),
      width = BiocGenerics::width(data),
      chromosome = paste0("chr", gsub('chr','', chr)),
      id = if (!isTRUE(tad.id)) {paste0(round(BiocGenerics::width(data) / 1e3), "Kb")} else {names(data)},
      genome = "custom",
      name = names(tad.lst)[i],
      stacking = "dense",
      featureAnnotation = "id",
      fontcolor.feature = "black",
      filled.contour = "black",
      col = "black",
      fill = colors.lst[i],
      cex = 0.6
    )
    tadTracks <- append(tadTracks, Track)
  }

  ##############################
  # bigwigTracks
  ##############################
  if (.Platform$OS.type == "windows") {bigwigPath.lst = NULL}
  if (is.null(bigwigPath.lst)) {
    bigwigTracks <- NULL
  } else {

    if (!is.list(bigwigPath.lst)) {
      stop("bigwigPath.lst must be a list")
    }

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

    bigwigTracks <- NULL

    # creation of the import region as GRange
    bin.gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges = IRanges::IRanges(start = start, end = stop),
      strand = "*"
    )
    GenomeInfoDb::seqlengths(bin.gr) <- GenomeInfoDb::seqlengths(tad.lst[[1]])[as.character(chr)]

    GenomeInfoDb::seqlevels(bin.gr) <- as.character(bigwig.chr)

    for (i in 1:length(bigwigPath.lst)) {

      # import bigwig files
      cov.gr <- rtracklayer::import.bw(bigwigPath.lst[[i]],
                                       which = bin.gr, as = "GRanges"
      )

      # Sanity check of chr names of bigwig file
      if (length(cov.gr[GenomeInfoDb::seqnames(cov.gr) == as.character(bigwig.chr)]) == 0) {
        stop(paste0("there is no chromosome ", bigwig.chr, " in the bigwig file"))
      }

      data.gr <- GenomicRanges::GRanges(
        seqnames = paste0("chr", gsub('chr','', chr)),
        ranges = IRanges::ranges(cov.gr),
        strand = "*",
        score = GenomicRanges::mcols(cov.gr)[1]
      )

      bigwigTrack <- Gviz::DataTrack(data.gr,
                                     chr = paste0("chr", gsub('chr','', chr)),
                                     type = "hist", aggregation = bigwig.xaxis,
                                     window = "fixed", windowSize = bigwig.binwidth,
                                     fill = colors.lst[i],
                                     name = names(bigwigPath.lst[i]),
                                     transformation = transformation.method
      )

      bigwigTracks <- append(bigwigTracks, bigwigTrack)
    }
  }

  #####################################
  # annotTracks
  #####################################
  if (is.null(annot.lst)) {
    annotTracks <- NULL
  } else {
    if (!is.list(annot.lst)) {
      stop("annot.lst must be a list")
    }

    annotTracks <- NULL
    for (i in 1:length(annot.lst)) {
      if (is.null(annot.col)) {
        data1 <- annot.lst[[i]][GenomeInfoDb::seqnames(annot.lst[[i]]) == chr &
                                  BiocGenerics::end(annot.lst[[i]]) >= start &
                                  GenomicRanges::start(annot.lst[[i]]) <= stop]

        if (length(data1) == 0) {
          annotTrack <- NULL
        } else {
          annotTrack <- Gviz::AnnotationTrack(
            start = GenomicRanges::start(data1),
            width = BiocGenerics::width(data1),
            chromosome = paste0("chr", gsub('chr','', chr)),
            strand = BiocGenerics::strand(data1),
            group = names(data1), groupAnnotation = "group",
            genome = "custom", name = names(annot.lst[i]), col = "deepskyblue4",
            col.line = "darkblue"
          )
        }
      } else {
        annotTrack <- NULL

        # Sanity check of annot.col
        if (length(GenomicRanges::mcols(annot.lst[[i]])) < annot.col) {
          stop(paste0("wrong annot.col number. There is only ", length(GenomicRanges::mcols(annot.lst[[i]])), " column(s) with metadata"))
        }

        data1 <- annot.lst[[i]][GenomeInfoDb::seqnames(annot.lst[[i]]) == chr &
                                  BiocGenerics::end(annot.lst[[i]]) >= start &
                                  GenomicRanges::start(annot.lst[[i]]) <= stop]

        if (length(data1) == 0) {
          annotTrack <- NULL
        } else {
          annotTrack <- Gviz::AnnotationTrack(
            start = GenomicRanges::start(data1),
            width = BiocGenerics::width(data1),
            chromosome = paste0("chr", gsub('chr','', chr)),
            strand = BiocGenerics::strand(data1),
            group = GenomicRanges::mcols(data1)[, annot.col], groupAnnotation = "group",
            genome = "custom", name = names(annot.lst[i]), col = "deepskyblue4",
            col.line = "darkblue", cex.feature = 0.5, cex.group = 0.7,
            just.group = "below"
          )
        }
      }
      annotTracks <- append(annotTracks, annotTrack)
    }
  }

  #####################################
  # bedgraphTracks
  #####################################

  #1 if bedgraphPath = NULL
  if (is.null(bedgraph.lst)) {
    bedgraphTracks <- NULL
  } else {

    #Stop if bedgraphPath is not a list
    if (!is.list(bedgraph.lst)) {
      stop("bedgraph.lst must be a list")
    }

    #If bedgraphPath[[1]] is not a list (list within list)
    if (!inherits(bedgraph.lst[[1]], "list")) {
      bedgraph.lst = list(bedgraph.lst)
      names(bedgraph.lst) = bedgraph.name
    }

    bedgraphTracks = NULL

    #For each tracks
    for (l in 1:length(bedgraph.lst)){
      data <- NULL

      #For each lines
      for (i in 1:length(bedgraph.lst[[l]])) {

        #read bedgraph datas
        ##if dataframe
        if (is.data.frame(bedgraph.lst[[l]][[i]])) {
          data0 = bedgraph.lst[[l]][[i]][,1:4]
        }
        ##if path
        if (is.character(bedgraph.lst[[l]][[i]])) {
          data0 = utils::read.table(bedgraph.lst[[l]][[i]], header = FALSE, sep = "\t")[,1:4]
        }
        ##if GRanges
        if (inherits(bedgraph.lst[[l]][[i]], "GRanges")) {
          data0 = as.data.frame(bedgraph.lst[[l]][[i]], row.names = NULL)[,c(1:3,6)]
        }

        #filter chr
        names(data0) = c(paste0("V", 1:4))
        data1 <- dplyr::filter(data0, V1 == chr)

        #If bedgraph values are all NA (ie no data in that chr)
        if (is.na(summary(data1$V4)[4])) {
          message(paste0(names(bedgraph.lst[l]), " datas of ", names(bedgraph.lst[[l]][i]), " is empty in that chr!"))
          next}

        #filter outliers
        if (bedgraph_outliers != 0) {
          data1 <- data1 %>% dplyr::filter(V4 < stats::quantile(V4, na.rm = TRUE, 1 - bedgraph_outliers)) %>%
            dplyr::filter(V4 > stats::quantile(V4, na.rm = TRUE, bedgraph_outliers))
        }

        #filter area (start & stop)
        data1 <- data1 %>% dplyr::filter(V3 >= start & V2 <= stop)

        if (is.na(summary(data1$V4)[4])) {
          message(paste0(names(bedgraph.lst[l])," datas of ", names(bedgraph.lst[[l]][i]), " is empty in that area!"))
          next}

        data1$V1 = as.character(data1$V1)
        names(data1) <- c("chr", "start", "end", names(bedgraph.lst[[l]][i]))
        data <- append(data, list(data1))
      }

      # if there is bedgraph datas in at least 1 line
      if (!is.null(data)) {

        # join dataframes, sort columns according to names (mcols names from GRange must be sorted) and create GRange
        data.gr <- data %>%
          Reduce(function(...) dplyr::full_join(..., by = c("chr", "start", "end")), .) %>%
          dplyr::select("chr", "start", "end", base::sapply(1:length(data), function(x){names(data[[x]][4])})) %>%
          dplyr::mutate(chr = paste0("chr", gsub('chr','', chr))) %>%
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

        bedgraphTracks.temp <- Gviz::DataTrack(
          range = sort(data.gr), na.rm = T,
          groups = names(GenomicRanges::mcols(data.gr)),
          type = "l", lwd = 2, name = names(bedgraph.lst[l]),
          col = colors.lst[1:length(bedgraph.lst[[l]])][order(names(GenomicRanges::mcols(data.gr)))]
        )

        bedgraphTracks = append(bedgraphTracks, bedgraphTracks.temp)
      }
    }
  }



  #####################################
  # grids
  #####################################
  data = GenomicRanges::restrict(sort(unlist(methods::as(tad.lst, "GRangesList"))),
                                 start = start, end = stop)
  borders = unique(sort(c(GenomicRanges::start(data[GenomeInfoDb::seqnames(data) == chr]),start, BiocGenerics::end(data[GenomeInfoDb::seqnames(data) == chr]),stop)))

  ht <- Gviz::HighlightTrack(trackList = c(areaTrack, tadTracks,
                                           bedgraphTracks, bigwigTracks,
                                           annotTracks),
                             start = borders,
                             end = borders,
                             chr = paste0("chr", gsub('chr','', chr)), col = "grey70", fill = "white")

  #####################################
  #Plot
  #####################################
  if (!is.null(title)) {main = title} else {main = ""}
  Gviz::plotTracks(c(ideoTrack, ht),
                   groupAnnotation = "group",
                   from = start,
                   to = stop,
                   background.title = "grey30", grid = TRUE, v = 0,
                   main = main
  )
}


