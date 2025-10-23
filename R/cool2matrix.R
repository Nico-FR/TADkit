#' @title Import matrix from cool or mcool file
#'
#' @description From cool or mcool files, `cool2matrix` import the interaction matrix for one chromosome as a `dgCMatrix` (upper triangular and sparse Matrix).
#' If balance = `TRUE`, `cool2matrix` returns the normalized counts.
#'
#' @details The cool file format is an efficient storage format for high resolution genomic interaction matrices, developed and maintained by the Mirny lab (https://github.com/open2c/cooler)
#' `cool2matrix` use the indexes provided by the cooler data model to extract the intra-chromosomal counts for a given chromosome.
#' As cool files store genomic interactions for only one resolution, `bin.width` must be set to `NA`.
#' While, as mcool files store genomic interactions for multiples resolutions, the chosen resolution must be set with the `bin.width` parameter.
#'
#' @param cool.path The full path of the cool or mcool file.
#' @param bin.width Bin width in base pair (i.e resolution) for mcool file only. Default is `NA` to parse cool file.
#' @param chr The selected chromosome.
#' @param balance Logical. Weather or not to use balanced counts instead of raw counts. Default = `FALSE`.
#' @param balancing_name Character that must correspond to the name given to the normalization method used. The most common names (for HiC normalisation) are "weight" (default), "KR", "VC".
#' @param verbose Logical. Whether or not to print messages. Default = `TRUE`.
#'
#' @return A `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @importFrom Matrix triu sparseMatrix
#' @importFrom rhdf5 h5read
#' @importFrom dplyr filter
#' @importFrom methods as
#' @importFrom magrittr %>%
#' @examples
#' # see vignette("Turorial_TADkit_R_package") or github (https://github.com/Nico-FR/TADkit)
#'
#' @export
#'
#'
cool2matrix <- function(cool.path, chr, bin.width = NA, balance = FALSE, balancing_name = "weight", verbose = TRUE) {

  if (verbose) {
    if (!is.na(bin.width)) {message("Parsing .mcool file.")} else {message("Parsing .cool file.")}
  }

  #mcool path
    uri <- function(cool.path) {
        if (!is.numeric(bin.width)) return(cool.path)
        return(
            paste(
                "resolutions",
                format(bin.width, scientific = FALSE),
                cool.path,
                sep = "/"
            )
        )
    }

    #verif bin.width
    if (!is.na(bin.width)) {

      if (!is.numeric(bin.width)) {
        bin.width = as.numeric(bin.width)
      }

      test_resolution = try(rhdf5::h5readAttributes(file = cool.path, name = uri("")), silent = TRUE)

      #if bin.width not available
      if (inherits(test_resolution, "try-error")) {

        #available resolutions
        ar = (rhdf5::h5ls(cool.path) %>% dplyr::filter(group == "/resolutions"))$name
        stop("\n '", bin.width, "' is not an available resolution.", " Resolutions available:\n", ar %>% as.numeric %>% sort %>% paste0(collapse = ", "), ".")
      }
    }

    # The list of available chromosomes
    chromosomes = rhdf5::h5read(file = cool.path, name = uri("chroms/name"))

    if (!(chr %in% chromosomes)) {
       stop("\n '", chr, "' is not a valid chromosome.", "\nChromosomes available: ", paste0(chromosomes, collapse = ", "), ".")
    }

    # Index of chrom in the chromosome list
    cid = match(chr, chromosomes)

    # chrom_offset: which row in the bin table each chromosome first appears (0-based)
    chrom_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/chrom_offset"))
    # chrom_lo: first occurence of chrom in the bin table 0-based
    chrom_lo = chrom_offset[cid] + 1   # 1-based
    # chrom_hi : last occurence of next chrom in the bin table 0-based)
    # chrom_offset[cid + 1] first occurence of next chrom (0-based)
    # in 1-based notation => last occurence of chrom in the bin table
    chrom_hi = chrom_offset[cid + 1]

    # bin1_offset:  which row in the pixel table each bin1 ID first appears.
    bin1_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/bin1_offset"))
    # lo first occurence of bin number chrom_lo in the pixel table (0-based)
    lo = bin1_offset[chrom_lo] + 1  # 1-based
    # hi first occurence of bin number chrom_hi + 1 in the pixels table (0-based)
    hi = bin1_offset[chrom_hi + 1]  # 1-based => last line with chrom_hi

    # chrom intra-chromosomal interactions
    id1 = rhdf5::h5read(file = cool.path, name = uri("pixels/bin1_id"), index = list(lo:hi))
    id2 = rhdf5::h5read(file = cool.path, name = uri("pixels/bin2_id"), index = list(lo:hi))
    interactions = rhdf5::h5read(file = cool.path, name = uri("pixels/count"), index = list(lo:hi))

    # Now we limit the interactions to intrachromosome interactions
    # hence removing bins starting from chrom_hi
    # chrom_hi corresponds to next chrom in 0-based
    # Let the indexes start at 0 (hence remove chrom_lo)
    i = id1[which(id2 < chrom_hi)] - chrom_lo + 1
    j = id2[which(id2 < chrom_hi)] - chrom_lo + 1
    x = interactions[which(id2 < chrom_hi)]

    #check if max(i) < nb bins (i.e is there a gap at the end of the chr?)
    #then, add the last bin indexes
    if (max(i) < chrom_hi - chrom_lo) {
      i = append(i, chrom_hi - chrom_lo)
      j = append(j, chrom_hi - chrom_lo)
      x = append(x, 0)
    }

    m = Matrix::sparseMatrix(i = i + 1, j = j + 1, x  =  as.numeric(x))

    if (balance) {

      # the list of available normalisation names
      tmp = as.data.frame(rhdf5::h5read(file = cool.path, name = uri("bins/"))) %>% dplyr::filter(chrom == chr) %>% names
      an = tmp[!tmp %in% c("start", "end", "chrom")]
      if (!(balancing_name %in% an)) {
        stop("\n '", balancing_name, "' normalisation is not available.", "\nNormalisations available are: ", paste0(an, collapse = ", "), ".")
      }

      if (verbose) {message("\nBalancing")}
      # Fetch the weights corresponding to the chromosome
      bins <- as.data.frame(rhdf5::h5read(file = cool.path, name = uri("bins/"))) %>% dplyr::filter(chrom == chr)
      w = bins[,balancing_name]
      #upper matrix weight for balancing
      mat_weight = Matrix::triu((w %*% t(w)))
      mat_weight[is.na(mat_weight)] <- 0 #remove NaN
      #cell by cell multiplication by the matrix weight
      mat = m * mat_weight
      if (!inherits(mat[1], "dgCMatrix")) {
        mat = methods::as(m * mat_weight, "CsparseMatrix")}
      return(mat)
    } else {
      return(m)
    }
}
