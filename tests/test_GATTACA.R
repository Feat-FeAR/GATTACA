source(file.path(ROOT, "src", "GATTACA.R"))

get_test_data("platinum_expr_set")
get_test_data("platinum_real_fcs")
get_test_data("previous_limma_out")

limma_out <- run_limma(
  platinum_expr_set, c(rep("A", 9), rep("B", 9)),
  contrasts = c("A-B", "B-A"),
  technical_replicates = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
  fc_threshold = 0.5
)

test_that("run_limma gives two dataframes", expect_length(limma_out, 2))

test_that("run_limma gives correct column names", {
  # Check colnames
  expect_setequal(
    colnames(limma_out$`A-B`),
    c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "markings")
  )
  expect_setequal(
    colnames(limma_out$`B-A`),
    c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "markings")
  )
})

test_that("run_limma preserves all genes in the output", {
  # Check rownames
  expect_setequal(
    rownames(limma_out$`A-B`), rownames(platinum_expr_set)
  )
  expect_setequal(
    rownames(limma_out$`B-A`), rownames(platinum_expr_set)
  )
})

test_that("run_limma gives the identical log fold changes", {
  # Check that the logFCs are unchanged from a previous run
  expect_equal(
    limma_out$`A-B`$logFC[order(rownames(limma_out$`A-B`))],
    previous_limma_out$`A-B`$logFC[order(rownames(limma_out$`A-B`))]
  )
  expect_equal(
    limma_out$`B-A`$logFC[order(rownames(limma_out$`B-A`))],
    previous_limma_out$`B-A`$logFC[order(rownames(limma_out$`B-A`))]
  )
})

test_that("run_limma respects the rules to assign markings", {
  # Check markings
  ## 1) The adj.p.val should be less than 0.05 for all marked genes
  expect_true(
    all(limma_out$`A-B`$adj.P.Val[limma_out$`A-B`$markings != 0] < 0.05)
  )
  expect_true(
    all(limma_out$`B-A`$adj.P.Val[limma_out$`B-A`$markings != 0] < 0.05)
  )

  ## 2) The abs FC should be greater than 0.5 for all marked genes
  expect_true(
    all(abs(limma_out$`A-B`$logFC[limma_out$`A-B`$markings != 0]) > 0.5)
  )
  expect_true(
    all(abs(limma_out$`B-A`$logFC[limma_out$`B-A`$markings != 0]) > 0.5)
  )

  ## 3) The logFC for genes marked as UP should be positive, and vice-versa
  expect_true(
    all(limma_out$`A-B`$logFC[limma_out$`A-B`$markings == 1] > 0)
  )
  expect_true(
    all(limma_out$`A-B`$logFC[limma_out$`A-B`$markings == -1] < 0)
  )
  expect_true(
    all(limma_out$`B-A`$logFC[limma_out$`B-A`$markings == 1] > 0)
  )
  expect_true(
    all(limma_out$`B-A`$logFC[limma_out$`B-A`$markings == -1] < 0)
  )
})

test_that("run_limma finds the correct up and downregulated spikes", {
  testable_limma <- lapply(
    limma_out, \(x){
      x[rownames(x) %in% intersect(rownames(x), platinum_real_fcs$probe_id),]
    }
  )

  # I make this into a log so that the comparisons are easier.
  platinum_real_fcs$real_fc_log <- platinum_real_fcs$real_FC
  # The logs break with zero, so i pull it to one, as it is the same as no
  # difference (log2 of 0) for our purposes
  platinum_real_fcs$real_fc_log[platinum_real_fcs$real_fc_log == 0] <- 1
  platinum_real_fcs$real_fc_log <- log2(platinum_real_fcs$real_fc_log)

  # This essentially tests that "A" is really the difference and "B" the baseline
  # The number of true positives should be greater than zero.
  expect_gt(
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    0
  )
  expect_gt(
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    0
  )
  ## The inverse should be true for the opposite run
  expect_gt(
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    0
  )
  expect_gt(
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    0
  )

  # We also expect that we get no real TPs if we do it the wrong way around.
  # NOTE:: This is not zero, but about zero, as some false positives can have
  # this error (even though it is really rare). So i put some small threshold.
  expect_lt(
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    20
  )
  expect_lt(
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    20
  )
  ## And the reverse...
  expect_lt(
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    20
  )
  expect_lt(
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    20
  )

  # Extra:: Print out the confusion matrix so that the tester can catch extra
  # possible errors
  conf_matrix <- matrix(nrow = 2, ncol = 2)
  rownames(conf_matrix) <- c("Real Positives", "Real Negatives")
  colnames(conf_matrix) <- c("Detected Positives", "Detected Negatives")
  # So, we have:
  # +--------+--------+
  # | 1,1 TP | 1,2 FN |
  # +--------+--------+
  # | 2,1 FP | 2,2 TN |
  # +--------+--------+

  cat("\nConfusion matrices for LIMMA")
  cat("\nThis is the confusion matrix for A vs B, considering B the baseline for UP genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings != 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings != 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log > 0]
    )
  print(conf_matrix)

  # Remake it for the downregulated case
  cat("\nThis is the confusion matrix for A vs B, considering B the baseline for DOWN genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings != -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_limma$`A-B`)[testable_limma$`A-B`$markings != -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )
  print(conf_matrix)

  # Remake it for B-A, so it's the inverse. It's identical to the other one,
  # but markings == 1 becomes == -1 and vice-versa (and so UP is DOWN)
  cat("\nThis is the confusion matrix for B - A, considering B the baseline for DOWN genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings != -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings != -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_FC > 1]
    )
  print(conf_matrix)

  # Remake it for the downregulated case
  cat("\nThis is the confusion matrix for B vs A, considering B the baseline for UP genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings != 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_limma$`B-A`)[testable_limma$`B-A`$markings != 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )
  print(conf_matrix)
})

###############################################################################

get_test_data("previous_rankprod_out")

{
  sink("/dev/null");
  rankprod_out <<- run_rankprod(
    platinum_expr_set, c(rep("A", 9), rep("B", 9)),
    contrasts = c("A-B", "B-A"),
    technical_replicates = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6),
    fc_threshold = 0.5
  )
  sink()
}


test_that("run_rankprod gives two dataframes", expect_length(rankprod_out, 2))

test_that("run_rankprod gives correct column names", {
  # Check colnames
  expect_setequal(
    colnames(rankprod_out$`A-B`),
    c(
      "gene.index", "logFC", "RP/Rsum.UP", "pfp.UP", "P.value.UP",
      "RP/Rsum.DOWN", "pfp.DOWN", "P.value.DOWN", "markings", "AveExpr"
      )
  )
  expect_setequal(
    colnames(rankprod_out$`B-A`),
    c(
      "gene.index", "logFC", "RP/Rsum.UP", "pfp.UP", "P.value.UP",
      "RP/Rsum.DOWN", "pfp.DOWN", "P.value.DOWN", "markings", "AveExpr"
    )
  )
})

test_that("run_rankprod preserves all genes in the output", {
  # Check rownames
  expect_setequal(
    rownames(rankprod_out$`A-B`), rownames(platinum_expr_set)
  )
  expect_setequal(
    rownames(rankprod_out$`B-A`), rownames(platinum_expr_set)
  )
})

test_that("run_rankprod gives the identical log fold changes", {
  # Check that the logFCs are unchanged from a previous run
  expect_equal(
    rankprod_out$`A-B`$logFC[order(rownames(rankprod_out$`A-B`))],
    previous_rankprod_out$`A-B`$logFC[order(rownames(rankprod_out$`A-B`))]
  )
  expect_equal(
    rankprod_out$`B-A`$logFC[order(rownames(rankprod_out$`B-A`))],
    previous_rankprod_out$`B-A`$logFC[order(rownames(rankprod_out$`B-A`))]
  )
})

test_that("run_rankprod respects the rules to assign markings", {
  # Check markings
  ## 1) The pfp should be less than 0.05 for all marked genes
  expect_true(
    all(rankprod_out$`A-B`$pfp.UP[rankprod_out$`A-B`$markings == 1] < 0.05)
  )
  expect_true(
    all(rankprod_out$`A-B`$pfp.DOWN[rankprod_out$`A-B`$markings == -1] < 0.05)
  )
  expect_true(
    all(rankprod_out$`B-A`$pfp.UP[rankprod_out$`B-A`$markings == 1] < 0.05)
  )
  expect_true(
    all(rankprod_out$`B-A`$pfp.DOWN[rankprod_out$`B-A`$markings == -1] < 0.05)
  )

  ## 2) The abs FC should be greater than 0.5 for all marked genes
  expect_true(
    all(abs(rankprod_out$`A-B`$logFC[rankprod_out$`A-B`$markings != 0]) > 0.5)
  )
  expect_true(
    all(abs(rankprod_out$`B-A`$logFC[rankprod_out$`B-A`$markings != 0]) > 0.5)
  )

  ## 3) The logFC for genes marked as UP should be positive, and vice-versa
  expect_true(
    all(rankprod_out$`A-B`$logFC[rankprod_out$`A-B`$markings == 1] > 0)
  )
  expect_true(
    all(rankprod_out$`A-B`$logFC[rankprod_out$`A-B`$markings == -1] < 0)
  )
  expect_true(
    all(rankprod_out$`B-A`$logFC[rankprod_out$`B-A`$markings == 1] > 0)
  )
  expect_true(
    all(rankprod_out$`B-A`$logFC[rankprod_out$`B-A`$markings == -1] < 0)
  )
})

test_that("run_rankprod finds the correct up and downregulated spikes", {
  testable_rankprod <- lapply(
    rankprod_out, \(x){
      x[rownames(x) %in% intersect(rownames(x), platinum_real_fcs$probe_id),]
    }
  )

  # I make this into a log so that the comparisons are easier.
  platinum_real_fcs$real_fc_log <- platinum_real_fcs$real_FC
  # The logs break with zero, so i pull it to one, as it is the same as no
  # difference (log2 of 0) for our purposes
  platinum_real_fcs$real_fc_log[platinum_real_fcs$real_fc_log == 0] <- 1
  platinum_real_fcs$real_fc_log <- log2(platinum_real_fcs$real_fc_log)

  # This essentially tests that "A" is really the difference and "B" the baseline
  # The number of true positives should be greater than zero.
  expect_gt(
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    0
  )
  expect_gt(
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    0
  )
  ## The inverse should be true for the opposite run
  expect_gt(
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    0
  )
  expect_gt(
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    0
  )

  # We also expect that we get no real TPs if we do it the wrong way around.
  # NOTE:: This is not zero, but about zero, as some false positives can have
  # this error (even though it is really rare). So i put some small threshold.
  expect_lt(
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    20
  )
  expect_lt(
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    20
  )
  ## And the reverse...
  expect_lt(
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    ),
    20
  )
  expect_lt(
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    ),
    20
  )

  # Extra:: Print out the confusion matrix so that the tester can catch extra
  # possible errors
  conf_matrix <- matrix(nrow = 2, ncol = 2)
  rownames(conf_matrix) <- c("Real Positives", "Real Negatives")
  colnames(conf_matrix) <- c("Detected Positives", "Detected Negatives")
  # So, we have:
  # +--------+--------+
  # | 1,1 TP | 1,2 FN |
  # +--------+--------+
  # | 2,1 FP | 2,2 TN |
  # +--------+--------+

  cat("\nConfusion matrices for RANKPROD")
  cat("\nThis is the confusion matrix for A vs B, considering B the baseline for UP genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings != 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log > 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings != 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log > 0]
    )
  print(conf_matrix)

  # Remake it for the downregulated case
  cat("\nThis is the confusion matrix for A vs B, considering B the baseline for DOWN genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings != -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings == -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_rankprod$`A-B`)[testable_rankprod$`A-B`$markings != -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )
  print(conf_matrix)

  # Remake it for B-A, so it's the inverse. It's identical to the other one,
  # but markings == 1 becomes == -1 and vice-versa (and so UP is DOWN)
  cat("\nThis is the confusion matrix for B - A, considering B the baseline for DOWN genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings != -1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_FC > 1]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings != -1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_FC > 1]
    )
  print(conf_matrix)

  # Remake it for the downregulated case
  cat("\nThis is the confusion matrix for B vs A, considering B the baseline for UP genes:\n")
  conf_matrix[1,1] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[1,2] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings != 1] %in%
        platinum_real_fcs$probe_id[platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,1] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings == 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )

  conf_matrix[2,2] <-
    sum(
      rownames(testable_rankprod$`B-A`)[testable_rankprod$`B-A`$markings != 1] %in%
        platinum_real_fcs$probe_id[! platinum_real_fcs$real_fc_log < 0]
    )
  print(conf_matrix)
})
