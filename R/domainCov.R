#' @title Distribution of annotation density within domains
#'
#' @description Distribution of annotations density (e.g. density  of genes) within domains (e.g. TADs, compartments...).
#' `domainCov()` calculates the annotations density per bin and measures the relative position of each bin within the domain
#'
#' @details As an example, this function take all TAD domains and calculate the density of annotations by:
#' * calculate the annotations density per bin,
#' * calculate the relative positions of each bin within his corresponding TAD,
#' * plot the smoothed bin density (or zscore bin density) according to relative positions.
#'
#' @param domain.gr `GRanges` with domains (TADs, compartments...).
#' @param annot.gr `GRanges` with genomic annotations (genes, repeat elements...).
#' @param annot.col Column number (metadata columns) used to split annotations. By default, annot.col = `NULL`. Use `"strand"` in order to separate annotations according to their strands.
#' @param domain.col Column number (metadata columns) used to split domain classes (e.g. column to differentiate compartments A and B). By default, domain.col is `NULL`.
#' @param bin.width Size of the bin in base pairs. It should match the bin.width of the matrix used to call the domains.
#' @param output Default is `"plot"` to return a `ggplot`. Use `"data"` to return the `dataframe` used to produce the plot.
#' @param norm Logical. Normalized density using relative content (zscore of density). Default is `FALSE`.
#'
#' @return `ggplot`
#'
#' @import GenomeInfoDb
#' @importFrom dplyr select filter
#' @import ggplot2
#' @import GenomicRanges
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @examples
#' # see vignette("Turorial_TADkit_R_package") or on github (https://github.com/Nico-FR/TADkit)
#'
#' library(EnsDb.Hsapiens.v86)
#' library(GenomicFeatures)
#'
#' # get gene annotations
#' genomic.gr =  genes(EnsDb.Hsapiens.v86, filter = ~ seq_name == c(1:22))
#' seqlevelsStyle(genomic.gr) = "UCSC"
#' genes.gr = genomic.gr[as.character(genomic.gr$gene_biotype) == "protein_coding"]
#'
#' # get compartments annotations
#' comp.gr = PC1calling(PC1_250kb.gr)
#'
#' domainCov(domain.gr = comp.gr, annot.gr = genes.gr,
#'           domain.col = 1,
#'           bin.width = 50e3, norm = TRUE)
#' @export
#'
domainCov <- function(domain.gr, annot.gr, domain.col = NULL, annot.col = NULL, bin.width = 50e3, norm = FALSE, output = "plot") {

  #annot.col parameter
  nb_metadatacolumns = length(GenomicRanges::mcols(annot.gr))
  if (is.null(annot.col)) {
    #mcols needed for analysis: add mcols with "all" for all annotations
    GenomicRanges::mcols(annot.gr)[,nb_metadatacolumns + 1] = "all"
    names(GenomicRanges::mcols(annot.gr))[nb_metadatacolumns + 1] = "type"
    annot.col = nb_metadatacolumns + 1
  }

  if (annot.col == "strand") {
    #add mcols with the stand of annotations
    GenomicRanges::mcols(annot.gr)[,nb_metadatacolumns + 1] = GenomicRanges::strand(annot.gr)
    names(GenomicRanges::mcols(annot.gr))[nb_metadatacolumns + 1] = "strand"
    annot.col = nb_metadatacolumns + 1
  }

  #sanity check of input parameters
  if (length(GenomicRanges::mcols(annot.gr)) < annot.col) {
    stop(paste0("Wrong annot.col number. There is only ", length(GenomicRanges::mcols(annot.gr)), " column(s) with metadatas."))
  }

  if (length(GenomicRanges::mcols(domain.gr)) < ifelse(is.null(domain.col), 0 ,domain.col)) {
    stop(paste0("Wrong domain.col number. There is only ", length(GenomicRanges::mcols(domain.gr)), " column(s) with metadatas."))
  }

  if (is.na(mean(GenomeInfoDb::seqlengths(domain.gr)))) {
    stop("There is no seqlenths in domain.gr, see dataframes2grange function")
  }

  #genome binning
  bin.gr = GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(domain.gr), tilewidth = bin.width, cut.last.tile.in.chrom = TRUE) %>%
    #keep domain.gr in chromosomes on annot.gr
    GenomeInfoDb::keepSeqlevels(seqlevels(annot.gr)[seqlevels(annot.gr) %in% seqlevels(domain.gr)],
                                pruning.mode = "coarse") %>% GenomeInfoDb::sortSeqlevels()

  #keep annot.gr in chromosomes on domain.gr
  annot2.gr = GenomeInfoDb::keepSeqlevels(annot.gr,
                                          seqlevels(domain.gr)[seqlevels(domain.gr) %in% seqlevels(annot.gr)],
                                          pruning.mode = "coarse") %>% GenomeInfoDb::sortSeqlevels()

  # add density for each annot.type (annot.col)
  annot.type = GenomicRanges::mcols(annot2.gr)[, annot.col] %>% as.character %>% unique
  for (c in annot.type) {
    bin.gr = GenomicRanges::binnedAverage(bins = bin.gr,
                                          numvar = GenomicRanges::coverage(annot2.gr[GenomicRanges::mcols(annot2.gr)[, annot.col] == c]),
                                          na.rm = FALSE, varname = c)
  }

  bin.gr <- GenomicRanges::resize(bin.gr, 1, fix = "center") #keep middle of each been

  #add domain names (ie row number) overlapped for each bin
  bin.gr$domainHit <- GenomicRanges::findOverlaps(bin.gr, domain.gr, select = "arbitrary")

  bin.gr = bin.gr[!is.na(bin.gr$domainHit)] #remove bin that do not overlap a domain

  # add domains datas needed (start, end, domain.col)
  bin.gr$domain_start = GenomicRanges::start(domain.gr[bin.gr$domainHit])
  bin.gr$domain_end = GenomicRanges::end(domain.gr[bin.gr$domainHit])
  if (is.numeric(domain.col)) {
    bin.gr$domain = GenomicRanges::mcols(domain.gr)[bin.gr$domainHit, domain.col]
  }

  #add relative position
  bin.gr$relative_position <- (GenomicRanges::start(bin.gr) - bin.gr$domain_start) / (bin.gr$domain_end - bin.gr$domain_start)

  #add zscore
  for (c in annot.type) {
    GenomicRanges::mcols(bin.gr)[, paste0("zscore_", c)] = scale(GenomicRanges::mcols(bin.gr)[, c])
  }

  if (output == "data") {
    return(bin.gr)
  }

  if (output == "plot") {

    data <- as.data.frame(bin.gr)
    names(data) = c("seqnames", "start", "end", "width", "strand", names(GenomicRanges::mcols(bin.gr)))

    if (is.null(domain.col)) {

      if (norm == FALSE) {
        data2 = data %>% dplyr::select(as.character(annot.type),"relative_position") %>% tidyr::gather(annot.col, "density", 1:length(annot.type))
        p1 = ggplot2::ggplot(data2, ggplot2::aes(y = density, x = relative_position, color = annot.col))+
          ggplot2::geom_smooth(se = FALSE)+
          ggplot2::labs(color = names(GenomicRanges::mcols(annot.gr[,annot.col])))
        return(p1)
      }

      if (norm == TRUE) {
        data4 <- data %>% dplyr::select(paste0("zscore_", annot.type), "relative_position")
        names(data4) <- c(as.character(annot.type), names(data4)[(length(annot.type)+1):length(data4)])
        data5 <- data4 %>% tidyr::gather(annot.col, "zscore", 1:length(annot.type))
        p1 = ggplot2::ggplot(data5, ggplot2::aes(y = zscore, x = relative_position, color = annot.col))+
          ggplot2::geom_smooth(se = FALSE)+
          ggplot2::labs(color = names(GenomicRanges::mcols(annot.gr[,annot.col])))
        return(p1)
      }

    }

    if (is.numeric(domain.col)) {

      if (norm == FALSE) {
        data2 = data %>% dplyr::select(as.character(annot.type), "relative_position","domain") %>% tidyr::gather(annot.col, "density", 1:length(annot.type))
        p1 = ggplot2::ggplot(data2, ggplot2::aes(y = density, x = relative_position, color = annot.col))+
          ggplot2::geom_smooth(se = FALSE)+ggplot2::facet_wrap(.~domain, scales = "free_y")+
          ggplot2::labs(color = names(GenomicRanges::mcols(domain.gr[, domain.col])))
        return(p1)
      }

      if (norm == TRUE) {
        data4 <- data %>% dplyr::select(paste0("zscore_", annot.type), "relative_position","domain")
        names(data4) <- c(as.character(annot.type), names(data4)[(length(annot.type)+1):length(data4)])
        data5 <- data4 %>% tidyr::gather(annot.col, "zscore", 1:length(annot.type))
        p3 = ggplot2::ggplot(data5, ggplot2::aes(y = zscore, x = relative_position, color = annot.col))+
          ggplot2::geom_smooth(se = FALSE)+ggplot2::facet_wrap(.~domain)+
          ggplot2::labs(color = names(GenomicRanges::mcols(annot.gr[,annot.col])))
        return(p3)
      }
    }
  }
}


























