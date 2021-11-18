source(file.path(ROOT, "src", "GATTACA.R"))

test_that("fix_replicates() returns the correct output", {
  input_data <- data.frame(
    "A.1" = c(1, 1, 1, 1),
    "A.2" = c(2, 3, 4, 5),
    "A.3" = c(10, 10, 10, 10),
    "A.4" = c(20, 30, 40, 50),
    "B.1" = c(100, 100, 100, 100),
    "B.2" = c(200, 300, 400, 500),
    "B.3" = c(1000, 1000, 1000, 1000),
    "B.4" = c(2000, 3000, 4000, 5000),
    row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
  )

  expect_equal(
    fix_replicates(
      input_data, groups = c("A", "A", "A", "A", "B", "B", "B", "B"),
      technical_replicates = c(1, 1, 2, 2, 3, 3, 4, 4)
    ),
    list(
      data = data.frame(
        "A.1" = c(1.5, 2, 2.5, 3),
        "A.2" = c(15, 20, 25, 30),
        "B.1" = c(150, 200, 250, 300),
        "B.2" = c(1500, 2000, 2500, 3000),
        row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
      ),
      groups = c("A", "A", "B", "B"),
      pairings = NULL
    )
  )

  expect_equal(
    fix_replicates(
      input_data, groups = c("P", "P", "P", "P", "K", "K", "K", "K"),
      technical_replicates = c(1, 1, 2, 2, 3, 3, 4, 4)
    ),
    list(
      data = data.frame(
        "P.1" = c(1.5, 2, 2.5, 3),
        "P.2" = c(15, 20, 25, 30),
        "K.1" = c(150, 200, 250, 300),
        "K.2" = c(1500, 2000, 2500, 3000),
        row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
      ),
      groups = c("P", "P", "K", "K"),
      pairings = NULL
    )
  )

  expect_equal(
    fix_replicates(
      input_data, groups = c("A", "A", "A", "A", "B", "B", "B", "B"),
      technical_replicates = c(1, 1, 2, 2, 3, 3, 4, 4),
      pairings = c(1,1,2,2,1,1,2,2)
    ),
    list(
      data = data.frame(
        "A.1" = c(1.5, 2, 2.5, 3),
        "A.2" = c(15, 20, 25, 30),
        "B.1" = c(150, 200, 250, 300),
        "B.2" = c(1500, 2000, 2500, 3000),
        row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
      ),
      groups = c("A", "A", "B", "B"),
      pairings = c(1, 2, 1, 2)
    )
  )
})

test_that("fix_replicates() crashes whith invalid inputs", {
  expect_error(
    fix_replicates(
      input_data, groups = c("P", "P", "K", "K"),
      technical_replicates = c(1, 1, 2, 2, 3, 3, 4, 4)
    )
  )
  expect_error(
    fix_replicates(
      input_data, groups = c("P", "P", "P", "P", "K", "K", "K", "K"),
      technical_replicates = c(1, 1, 2, 2, 3, 4, 4)
    )
  )
  expect_error(
    fix_replicates(
      input_data, groups = c("P", "P", "P", "P", "K", "K", "K", "K"),
      technical_replicates = c(1, 1, 2, 2, 3, 4, 4),
      pairings = c(1, 2)
    )
  )
  expect_error(
    fix_replicates(
      input_data, groups = c("P", "P", "P", "P", "K", "K", "K", "K"),
      technical_replicates = c(1, 1, 2, 2, 3, 4, 4),
      pairings = c(1,1,1,2,1,1,2,2)
    )
  )
})

test_that("fix_replicates() leaves single replicates untouched", {
  input_data <- data.frame(
    "A.1" = c(1, 1, 1, 1),
    "A.2" = c(10, 10, 10, 10),
    "A.3" = c(20, 30, 40, 50),
    "B.1" = c(100, 100, 100, 100),
    "B.2" = c(1000, 1000, 1000, 1000),
    "B.3" = c(2000, 3000, 4000, 5000),
    row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
  )

  expect_equal(
    fix_replicates(
      input_data, groups = c("A", "A", "A", "B", "B", "B"),
      technical_replicates = c(1, 2, 2, 3, 4, 4)
    ),
    list(
      data = data.frame(
        "A.1" = c(1, 1, 1, 1),
        "A.2" = c(15, 20, 25, 30),
        "B.1" = c(100, 100, 100, 100),
        "B.2" = c(1500, 2000, 2500, 3000),
        row.names = c("gene_a", "gene_b", "gene_c", "gene_d")
      ),
      groups = c("A", "A", "B", "B"),
      pairings = NULL
    )
  )
})


test_that("run_limma gives expected output", {
  get_test_data("platinum_expr_set")
  get_test_data("platinum_real_fcs")

  limma_out <- run_limma(
    platinum_expr_set, c(rep("A", 9), rep("B", 9)),
    contrasts = c("A-B", "B-A"),
    technical_replicates = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6)
  )



})