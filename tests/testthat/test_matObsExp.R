library(testthat)
library(TADkit)
library(Matrix)

test_that("matObsExp basic functionality and optimization", {
  # Load the example dataset from TADkit
  env = new.env()
  load(system.file("data", "mat_HCT116_chr19_50kb.rda", package = "TADkit"), envir = env)
  # If system.file didn't find it (e.g. package not installed yet), load from local data/
  if (length(ls(env)) == 0) {
    load("../../data/mat_HCT116_chr19_50kb.rda", envir = env)
  }
  mat = env[["mat_HCT116_chr19_50kb"]]

  # Test that default output is OE and returns a CsparseMatrix (or dgCMatrix)
  res_oe = matObsExp(mat, "OE")
  expect_true(inherits(res_oe, "CsparseMatrix"))
  expect_equal(dim(res_oe), dim(mat))

  # Test that E output returns a dense matrix (Toeplitz)
  res_e = matObsExp(mat, "E")
  expect_true(inherits(res_e, "matrix") || inherits(res_e, "Matrix"))
  expect_equal(dim(res_e), dim(mat))

  # Test with a standard matrix input
  dense_mat = as.matrix(mat[1:100, 1:100])
  res_dense_oe = matObsExp(dense_mat, "OE")
  expect_true(inherits(res_dense_oe, "CsparseMatrix"))

  # Verify correctness of calculations by checking some specific non-zero elements
  # Find some non-zero elements
  s = Matrix::summary(mat)
  # Ensure we have at least one element with dist > 0
  non_diag = s[s$i != s$j, ]
  if (nrow(non_diag) > 0) {
    row_idx = non_diag$i[1]
    col_idx = non_diag$j[1]
    val = non_diag$x[1]
    dist_val = abs(row_idx - col_idx)

    # Expected value for this distance: average of all non-zero values at this distance
    all_non_zeros_at_dist = s[abs(s$i - s$j) == dist_val & s$i <= s$j, "x"]
    expected_mean = mean(all_non_zeros_at_dist)

    # OE ratio should be val / expected_mean
    expect_equal(res_oe[row_idx, col_idx], val / expected_mean)
  }
})

test_that("matObsExp works correctly and fast on large mock sparse matrices", {
  # Create a mock sparse matrix of large dimensions (e.g., 2000 x 2000)
  # which would have caused issues if grid expansion of 4 million cells was slow or heavy.
  # This serves as a scalability and performance safeguard.
  n = 2000
  row_indices = c(1, 10, 50, 100, 500, 1000, 1500, 1999)
  col_indices = c(2, 20, 100, 200, 600, 1100, 1600, 2000)
  values = c(1.5, 2.3, 4.1, 0.8, 12.0, 5.5, 3.2, 9.1)

  large_mat = Matrix::sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = values,
    dims = c(n, n)
  )

  # Should run instantly and successfully without any memory errors
  start_time = Sys.time()
  res_large = matObsExp(large_mat, "OE")
  end_time = Sys.time()

  expect_true(inherits(res_large, "CsparseMatrix"))
  expect_equal(dim(res_large), c(n, n))
  expect_true(as.numeric(end_time - start_time, units = "secs") < 2.0)
})
