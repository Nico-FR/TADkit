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
cool2matrix <- function(cool.path, chr, bin.width = NA, balance = FALSE, balancing_name = "weight", verbose = TRUE, chr2 = NULL) {

  ###################
  #cool.path = "~/mnt/genome3D/Var_struc/Mouton_onCow/Contact_maps/Distiller_Arima/norm/mcool/mouton8045.ARS-UCD1.2.mapq_10.4000.mcool"
  #chr = 1
  #chr2 = NULL
  #bin.width = 64000
  #balance = TRUE; balancing_name = "weight"; verbose = TRUE
  ###################

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

      #make sure bin.width is numeric
      if (!is.numeric(bin.width)) {
        bin.width = as.numeric(bin.width)
      }

      test_resolution = try(rhdf5::h5readAttributes(file = cool.path, name = uri("")), silent = TRUE) # test if the resolution exists

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
    if (!is.null(chr2)) {
      if (!(chr2 %in% chromosomes)) {
         stop("\n '", chr2, "' is not a valid chromosome.", "\nChromosomes available: ", paste0(chromosomes, collapse = ", "), ".")
      }
    }

    rank_1 = match(chr, chromosomes) #rank of chr in chromosomes list
    if (!is.null(chr2)) {
      rank_2 = match(chr2, chromosomes) #rank of chr2 in chromosomes list

      #check that rank_1 < ran_2
      if (!rank_1 < rank_2) { #need to switch chr order to match cool file
        flip.mat = TRUE #indicate that the matrix need to be flipped at the end
        rank_1.tmp = rank_1
        rank_1 = rank_2
        rank_2 = rank_1.tmp
      } else {
        flip.mat = FALSE
      }
    }

    # Fetch first bin IDs of chr
    chrom1_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/chrom_offset"))[rank_1:(rank_1 + 1)] + 1 #first bin ID (ie bin number) of chr & first bin ID of next chromosome
    chrom1_first_bin_index = rhdf5::h5read(file = cool.path, name = uri("indexes/bin1_offset"))[chrom1_offset[1]] + 1 #first occurrence (ie line number = index) of first bin of chr in pixel table

    # Fetch last bin IDs of chr or chr2
    if (is.null(chr2)) { #if only one chromosome
      chrom2_offset = chrom1_offset
      chrom2_last_bin_index = rhdf5::h5read(file = cool.path, name = uri("indexes/bin1_offset"))[chrom1_offset[2] - 1] #first occurrence (ie line number = index) of last bin of chr in pixel table
    } else { #if two chromosomes
      chrom2_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/chrom_offset"))[rank_2:(rank_2 + 1)] + 1 #first bin ID (ie bin number) of chr2 & first bin ID of next chromosome
      chrom2_last_bin_index = rhdf5::h5read(file = cool.path, name = uri("indexes/bin1_offset"))[chrom2_offset[2] - 1] #first occurrence (ie line number = index) of last bin of chr2 in pixel table
    }

    # number of bins for each chromosome
    nb_bin_chr1 = chrom1_offset[2] - chrom1_offset[1]
    nb_bin_chr2 = chrom2_offset[2] - chrom2_offset[1]

    if (verbose) {
      size_message = if(flip.mat) {paste0(nb_bin_chr2, "x", nb_bin_chr1)} else {paste0(nb_bin_chr1, "x", nb_bin_chr2)}
      if (!is.na(bin.width)) {
        message(
          paste0("Parsing .mcool file as ", size_message, " matrix at ", bin.width, " bp resolution."))}
      else {message(paste0("Parsing .cool file as ", size_message, " matrix at ", bin.width, " bp resolution."))}
    }

    # Read the pixels data (3 col tables) for chr
    bin1_id = rhdf5::h5read(file = cool.path, name = uri("pixels/bin1_id"), index = list(chrom1_first_bin_index:chrom2_last_bin_index))
    bin2_id = rhdf5::h5read(file = cool.path, name = uri("pixels/bin2_id"), index = list(chrom1_first_bin_index:chrom2_last_bin_index))
    interactions = rhdf5::h5read(file = cool.path, name = uri("pixels/count"), index = list(chrom1_first_bin_index:chrom2_last_bin_index))

    melted.mat = data.frame(bin1_id = bin1_id + 1,
                            bin2_id = bin2_id + 1, #convert to 1based bin IDs
                            count = interactions) %>%
      dplyr::filter(bin1_id %in% c(chrom1_offset[1]:(chrom1_offset[2] - 1))) %>% # filter only chr bins
      dplyr::filter(bin2_id %in% c(chrom2_offset[1]:(chrom2_offset[2] - 1))) # filter only chr2 bins

    melted.mat = melted.mat %>% dplyr::mutate(bin1_id = bin1_id - chrom1_offset[1] + 1, bin2_id = bin2_id - chrom2_offset[1] + 1) #reindex bins IDs to start at 1 for each chromosome

    #check if max(melted.mat$bin1_id) < nb bins (i.e is there a gap at the end of the chr?)
    #then, add the last bin index
    if (max(melted.mat$bin1_id) < nb_bin_chr1 | max(melted.mat$bin2_id) < nb_bin_chr2) {
      melted.mat = dplyr::bind_rows(
        melted.mat,
        data.frame(
          bin1_id = nb_bin_chr1,
          bin2_id = nb_bin_chr2,
          count = 0
        )
      )
    }

    # flip the matrix if needed
    if (!is.null(chr2) & flip.mat) {
      melted.mat = melted.mat %>%
        dplyr::rename(bin1_id = bin2_id,
                      bin2_id = bin1_id)
    }

    # Create sparse matrix
    m = Matrix::sparseMatrix(i = melted.mat$bin1_id,
                             j = melted.mat$bin2_id,
                             x = as.numeric(melted.mat$count))

    if (balance) {

      if (verbose) {message("\nBalancing")}

      # check normalization names (balancing_name)
      test_balancing_name = try(rhdf5::h5readAttributes(file = cool.path, name = uri(paste0("bins/", balancing_name))), silent = TRUE)
      if (inherits(test_balancing_name, "try-error")) {
        tmp = rhdf5::h5read(file = cool.path, name = uri("bins")) %>% names
        an = tmp[!tmp %in% c("start", "end", "chrom")]
        stop("\n '", balancing_name, "' normalisation is not available.", "\nNormalisations available are: ", paste0(an, collapse = ", "), ".")
        }

      # Fetch the normalization values corresponding to the chr
      w_chr1 = rhdf5::h5read(file = cool.path, name = uri("bins/weight"), index = list(chrom1_offset[1]:(chrom1_offset[2] - 1)))

      #normalization matrix
      if (is.null(chr2)) { #if only one chromosome
        #upper matrix weight for intra chromosomal balancing
        mat_weight = Matrix::triu(w_chr1 %*% t(w_chr1))
        } else { #if two chromosomes
          w_chr2 = rhdf5::h5read(file = cool.path, name = uri("bins/weight"), index = list(chrom2_offset[1]:(chrom2_offset[2] - 1)))
          #full matrix weight for inter chromosomal balancing
          mat_weight = Matrix::Matrix(w_chr1 %*% t(w_chr2))

          if (flip.mat) {
            mat_weight = t(mat_weight)
          }
        }

      mat_weight[is.na(mat_weight)] <- 0 #remove NaN

      #Balancing: cell by cell multiplication by the matrix normalization
      mat = m * mat_weight

      if (!inherits(mat, "dgCMatrix")) {
        mat = methods::as(mat, "CsparseMatrix")}
      return(mat)
    }

    return(m)
}
