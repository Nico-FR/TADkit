#' @title Orientation of PC1 score
#'
#' @description Assuming that A compartments have a higher expression than B compartments.
#' This function allows the PC1 score to be oriented from gene expression data.
#' Then the compartments A (actif) and B (inactif) can be called accurately.
#'
#' @details This function is divided into 3 steps:
#' * Compartments calling from the score (see PC1calling function).
#' * Measures median gene expression in compartments A and B per chromosome.
#' * Reverse the score if median expression of A is smaller than B compartments.
#'
#' Return a list with 3 `dataframe`.
#' * bedgraph as `GRanges` object with the score oriented.
#' * `dataframe` with expression of each gene (before orientation).
#' * `dataframe` with median expression of each compartment (before orientation).
#'
#' @param bedgraph.gr `GRanges` file with the score to be used for compartment calling (i.e PC1 values).
#' @param annot.gr `GRanges` file with gene annotations.
#' @param expression.data.frame `dataframe` with 3 columns:
#' * 1: gene IDs: must be identical IDs than `names(annot.gr)`,
#' * 2: common name for genes: could be the gene ID,
#' * 3: expression of each gene: such as raw count or `log(raw_count + 1)`.
#'
#' @return S3 class object with 3 `dataframes`.
#'
#' @import GenomeInfoDb
#' @importFrom dplyr select group_by summarise arrange filter rename
#' @importFrom tidyr spread
#' @importFrom stats median
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' # see vignette("Turorial_TADkit_R_package") or on github (https://github.com/Nico-FR/TADkit)
#'
compOrientation <- function(bedgraph.gr, annot.gr, expression.data.frame) {

  #local variables:
  comp <- chr <- med_exp <- A <- B <- . <- NULL

  #sanity check
  if (is.null(names(annot.gr))) {
    stop("annot.gr object must have names (i.e gene IDs)")
  }
  if (length(annot.gr[duplicated(names(annot.gr))]) > 0) {
    stop("annot.gr object must have unique names (i.e gene IDs)")
  }
  if (length(expression.data.frame[duplicated(expression.data.frame[,1]),1]) > 0) {
    stop("expression.data.frame[,1] must have unique IDs (i.e gene IDs)")
  }
  if (length(expression.data.frame) < 3) {
    stop("expression.data.frame must have 3 columns")
  }
  if (length(expression.data.frame) > 3) {
    warning("expression.data.frame have more than 3 columns, only the third column is used\n")
  }
  if (!inherits(bedgraph.gr, "GRanges")) {
    stop("bedgraph.gr must be a GRange object")
  }
  if (is.na(mean(seqlengths(bedgraph.gr), na.rm=T))) {
    stop("bedgraph.gr must have seqlengths datas (see dataframes2grange function)")
  }
  if (length(S4Vectors::mcols(bedgraph.gr)) != 1) {
    warning("bedgraph.gr object should have 1 metadata column, additional columns are ignored!")
  }

  #call compartments
  comp.gr = suppressWarnings(PC1calling(bedgraph.gr))

  #keep comp.gr in chromosomes from annot.gr
  undirected_comp.gr = comp.gr %>% GenomeInfoDb::keepSeqlevels(seqlevels(annot.gr)[seqlevels(annot.gr) %in% seqlevels(bedgraph.gr)],
                                                               pruning.mode = "coarse") %>% GenomeInfoDb::sortSeqlevels()

  #keep expression.data.frame IDs from annot.gr IDs
  expression.data.frame2 = expression.data.frame[expression.data.frame[,1] %in% names(annot.gr),1:3]

  #keep annot.gr in chromosomes from comp.gr
  annot.gr = GenomeInfoDb::keepSeqlevels(annot.gr,
                                         seqlevels(comp.gr)[seqlevels(comp.gr) %in% seqlevels(annot.gr)],
                                         pruning.mode = "coarse") %>% GenomeInfoDb::sortSeqlevels()

  # add unique names to each domain
  names(undirected_comp.gr) = 1:length(undirected_comp.gr)

  # position of the genes in relation to each compartment
  geneTopo <- suppressWarnings(
    geneTADtopo(domain.gr = undirected_comp.gr, annot.gr = annot.gr, expression.data.frame = expression.data.frame2, ifoverlap = "remove")
  )

  # merge gene, domain and expression datas
  undirected_comp.df = undirected_comp.gr %>% as.data.frame() %>% dplyr::select(seqnames, comp)
  undirected_comp.df$TAD_id = rownames(undirected_comp.df)
  rownames(undirected_comp.df) = NULL

  data = merge(geneTopo$geneTADtopo, undirected_comp.df)

  # add comp and chr to expression.df
  for (i in 1:length(data[,1])) {
    if (nrow(data$gene_exp.lst[[i]]) == 0) {next} #ignore genes without expression datas
    data$gene_exp.lst[[i]]$comp = data$comp[i]
    data$gene_exp.lst[[i]]$chr = data$seqnames[i]
  }

  # merge expression.df list
  data = do.call(rbind, data$gene_exp.lst) %>%
    dplyr::rename("exp" = 3)

  # stats
  medianExp = data %>% dplyr::group_by(comp, chr) %>% dplyr::summarise(med_exp = median(exp)) %>% dplyr::arrange(chr) %>% tidyr::spread(., comp, med_exp) %>% as.data.frame()
  medianExp$to_be_inverted = ifelse(medianExp$A < medianExp$B, -1, 1)

  message(paste0(length(medianExp[medianExp$to_be_inverted == -1,1]), "/", length(medianExp[,1]), " chromosomes have been reversed"))

  # if exp A == exp B ==> warning
  if (length((medianExp %>% filter(A == B))[,1])) {
    message(paste0("median expression levels of chromosome ", ((medianExp %>% filter(A == B))[,1]), " are identical between compartments A and B\n"))
  }

  #reverse score
  bedgraph_oriented = NULL
  for (c in unique(medianExp$chr)) {
    tmp = as.data.frame(bedgraph.gr[GenomeInfoDb::seqnames(bedgraph.gr) == c])
    tmp[,6] = tmp[,6] * medianExp$to_be_inverted[medianExp$chr == c]
    bedgraph_oriented = rbind(bedgraph_oriented, tmp)
  }
  bedgraph_oriented.gr = GenomicRanges::makeGRangesFromDataFrame(bedgraph_oriented, keep.extra.columns = T)

  # add chromosomes lengths...
  suppressWarnings(
    seqlengths(bedgraph_oriented.gr)[sort(names(seqlengths(bedgraph_oriented.gr)))] <- seqlengths(bedgraph.gr)[sort(names(seqlengths(bedgraph.gr)))])

  output = list()
  output$bedgraph_oriented = bedgraph_oriented.gr
  output$expression = data %>% `colnames<-`(c(names(data)[1:2], names(expression.data.frame2)[3], names(data)[4:5]))
  output$medianExp = medianExp[,1:3]
  return(output)
}








