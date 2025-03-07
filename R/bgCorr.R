#' @title Correlation between bedgraph scores
#'
#' @description Computes the correlations between bin scores (e.g. between insulation scores...) from 2 or more individuals.
#'
#'
#' @param bedgraph.lst list of `data.frame`, `GRanges` or full path of the bedgraph files (data frame without header and 4 columns tab separated) containing a score for each bin.
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param rm_chr chromosome to filter, default = "X".
#' @param Qnorm perform quantile normalisation beetween scores, default = FALSE.
#'
#' @return matrix with correlation values
#'
#' @import GenomicRanges
#' @importFrom utils read.table
#' @importFrom dplyr select full_join filter mutate_at
#' @importFrom preprocessCore normalize.quantiles
#' @examples
#' bedgraphs.lst = list(ind1 = IS_HCT116_chr19_5kb.bedgraph, ind2 = IS_HCT116_chr19_5kb.bedgraph)
#' bgCorr(bedgraphs.lst)
#'
#' @export
#'
bgCorr <- function(bedgraph.lst, method = "pearson", rm_chr = "X", Qnorm = FALSE) {

  . <- NULL

  #Stop if bedgraphPath is not a list
  if (!is.list(bedgraph.lst)) {
    stop("bedgraph.lst must be a list")
  }

  #read bedgraph datas
  ##if dataframe
  if (is.data.frame(bedgraph.lst[[1]])) {
    data1 =  base::lapply(bedgraph.lst, function(bg){
      bg[,1:4] %>% mutate_at(4, as.numeric) %>% mutate_at(1, as.character)})
  }
  ##if path
  if (is.character(bedgraph.lst[[1]])) {
    data1 = base::lapply(bedgraph.lst, function(bg){
      utils::read.table(bg, header = FALSE, sep = "\t")[,1:4] %>% mutate_at(4, as.numeric) %>% mutate_at(1, as.character)})
  }
  ##if GRanges
  if (inherits(bedgraph.lst[[1]], "GRanges")) {
    data1 = base::lapply(bedgraph.lst, function(bg){
      as.data.frame(bg) %>% dplyr::select(c(1:3,6)) %>% mutate_at(4, as.numeric)})
  }

  ##if names is null
  if (is.null(names(bedgraph.lst))) {
    names(data1) = 1:length(names(bedgraph.lst))
  }

  #add name to V4
  for (c in names(data1)) {
    names(data1[[c]]) = c("seqnames", "start", "end", c)
  }

  #merge bedgraphs
  data2 = data1 %>% Reduce(function(...) dplyr::full_join(..., by = c("seqnames", "start", "end")), .) %>% dplyr::filter(!seqnames %in% rm_chr) %>% select(-c("seqnames", "start", "end")) %>% as.matrix

  #quantile normalisation
  if (Qnorm == TRUE) {
    preprocessCore::normalize.quantiles(data2, copy=FALSE, keep.names=TRUE)
  }

  #create matrix
  mat=matrix(nrow = length(bedgraph.lst), ncol = length(bedgraph.lst))
  colnames(mat) = rownames(mat) = names(bedgraph.lst)

  for (ind1 in names(bedgraph.lst)) {
    for (ind2 in names(bedgraph.lst)) {
      if (ind1 == ind2) {next}
      mat[ind1, ind2] = stats::cor(data2[, ind1], data2[, ind2], use = "complete.obs", method = method)
    }
  }
  diag(mat) = 1

  return(mat)
}
