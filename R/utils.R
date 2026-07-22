# Silence "no visible binding for global variable" warnings from R CMD check
utils::globalVariables(c(
  "i", "e2", "s2", "start1", "end1", "start2", "end2", "chr1", "chr2", "x", "y", "e", "s", "chr", # MATplot
  "i", "e2", "s2", "start1", "end1", "start2", "end2", "chr1", "chr2", "x", "y", "e", "s", "chr", # mMATplot
  "j", # MATfeatures
  "density", "zscore", # areaCov
  ".", # bgCorr
  "comp", ".", "med_exp", "A", "B", # compOrientation
  "density", "relative_position", "zscore", # domainCov
  "mixed_position", # domainHist
  "nameHit", "nb_genes", "gene_strands", "nb_TADs", # geneTADtopo
  "j", ".", "expected", # matExpDist : <anonymous>
  "expected", # matExpDist
  "j", ".", "dist", # matObsExp
  ".", # viewPointInteract : <anonymous>
  ".", "mean_interact", # viewPointInteract
  'tad.chr', 's_y', 'e_y', 's_x', 'e_x', 'c_y', 'c_x', #interMATplot
  "bin1_id", "bin2_id"
))
