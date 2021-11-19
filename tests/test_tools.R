test_that("str_intersection() gives expected output", {
  expect_equal(
    str_intersection("abbac", "abd"), c("a", "b")
  )
  expect_equal(
    str_intersection("Aa", "a"), c("a")
  )
  expect_equal(
    str_intersection("", "abdcde"), c()
  )
  expect_equal(
    str_intersection("123))0", "12)0002"), c("1", "2", ")", "0")
  )
  expect_equal(
    str_intersection("cba", "abc"), c("c", "b", "a")
  )
})


test_that("find_nums() gives expected output", {
  expect_equal(
    find_nums("aandc12,3o0.2//23068"), c(12, 3, 0, 2, 23068)
  )
  expect_equal(
    find_nums(""), c()
  )
  expect_equal(
    find_nums(12330002), 12330002
  )
  expect_equal(
    find_nums("001,002"), c(1, 2)
  )
})


test_that("subtract.str() gives expected output", {
  expect_equal(
    subtract.str("ABab1", "ABBA is my favourite band!!11!."),
    " is my fvourite nd!!!."
  )
  expect_equal(
    subtract.str("", "ABBA is my favourite band!!11!."),
    "ABBA is my favourite band!!11!."
  )
  expect_equal(
    subtract.str(" ", "ABBA is my favourite band!!11!."),
    "ABBAismyfavouriteband!!11!."
  )
  expect_equal(
    subtract.str("ABBA is my favourite band!!11!.", ""),
    ""
  )
})


get_test_data("test_expression_data")
get_test_data("qq_norm_expr_data")

test_that("qq_normalize() gives expected output", {
  expect_equal(
    qq_normalize(test_expression_data), qq_norm_expr_data
  )
})


test_that("make.unique_from1() gives expected output", {
  expect_equal(
    make.unique_from1(
      c("banana", "banana", "papaya", "guava", "banana", "guava")
    ),
    c("banana.1", "banana.2", "papaya.1", "guava.1", "banana.3", "guava.2")
  )
  expect_equal(
    # This is retarded, but `make.unique` does this too.
    make.unique_from1(""), ".1"
  )
  expect_equal(
    make.unique_from1(c(1, 2, 3)), c("1.1", "2.1", "3.1")
  )
})


test_that("get_median() gives expected output", {
  expect_equal(
    get_median(data.frame(
      sample.1 = c(1, 2, 3),
      sample.2 = c(2, 2, 1),
      sample.3 = c(3, 4, 2),
      row.names = c("a", "b", "c")
    )),
    data.frame(
      Median = c(2, 2, 2),
      row.names = c("a", "b", "c")
    )
  )
})


test_that("bin_mean() gives expected output and errors", {
  expect_equal(
    bin_mean(c(1,2,3,4,5,6), n_bins = 3),
    c(1.5, 3.5, 5.5)
  )
  expect_equal(
    bin_mean(c(1,2,3,4,5,6,7), n_bins = 3),
    c(2, 4, 6)
  )
  expect_equal(
    bin_mean(1, n_bins = 1), 1
  )
  expect_equal(bin_mean(1:1000, n_bins = 1), mean(1:1000))
  expect_error(bin_mean(c()))
  expect_error(bin_mean(c(1,2,3), n_bins = 10))
  expect_error(bin_mean(c(1,2,3), bin_by = c(1, 2, 3)))
})


test_that("find_BH_critical_p() gives expected output and errors", {
  expect_equal(
    find_BH_critical_p(c(0.05, 0.01, 0.25)),
    0.03333333
  )
  expect_equal(
    find_BH_critical_p(c(0.25, 1, 0.55)),
    0.01666667
  )
  expect_equal(
    find_BH_critical_p(c(0.01, 0.01, 0.01)),
    0.06666667
  )
  expect_equal(
    find_BH_critical_p(c(0.01, 0.01, 0.01), 0),
    0
  )
  expect_error(find_BH_critical_p(c()))
  expect_error(find_BH_critical_p(c(1, 0.01, 0.12), -1))
  expect_warning(find_BH_critical_p(c(-1, 2, 0)))
})


# Test the design parser

test_that("Expansion Functions expand correctly", {
  expect_equal(
    expand_intercalated("banana,papaya", 2, 1),
    c("banana,papaya,banana,papaya")
  )
  expect_equal(
    expand_intercalated("banana,papaya", 1, 1),
    c("banana,papaya")
  )
  expect_equal(
    expand_intercalated("banana", 3, 1),
    c("banana,banana,banana")
  )
  expect_equal(
    expand_intercalated("banana*", 3, 1),
    c("banana1,banana2,banana3")
  )
  expect_equal(
    expand_intercalated("banana*,papaya*", 2, 3),
    c("banana3,papaya3,banana4,papaya4")
  )
  expect_error(
    expand_intercalated("", 10, 1)
  )
  expect_error(
    expand_intercalated(",banana", 2, 1)
  )
  expect_error(
    expand_intercalated("banana,,papaya", 2, 1)
  )

  expect_equal(
    expand_ordered("banana,papaya", 2, 1),
    c("banana,banana,papaya,papaya")
  )
  expect_equal(
    expand_ordered("banana,papaya", 1, 1),
    c("banana,papaya")
  )
  expect_equal(
    expand_ordered("banana", 3, 1),
    c("banana,banana,banana")
  )
  expect_equal(
    expand_ordered("banana*", 3, 1),
    c("banana1,banana2,banana3")
  )
  expect_equal(
    expand_ordered("banana*,papaya*", 2, 3),
    c("banana3,banana4,papaya3,papaya4")
  )
  expect_error(
    expand_ordered("", 10, 1)
  )
  expect_error(
    expand_ordered(",banana", 2, 1)
  )
  expect_error(
    expand_ordered("banana,,papaya", 2, 1)
  )
})

test_that("get_captures() gives expected output", {
  expect_equal(
    get_captures(
      "([A-Z])\\w+",
      "RegExr was created by gskinner.com, and is proudly hosted by Media Temple."
    ),
    c("RegExr", "R")
  )
  expect_equal(
    get_captures(
      "(RegExr)+",
      "RegExr was created by gskinner.com, and is proudly hosted by Media Temple."
    ),
    c("RegExr", "RegExr")
  )
  expect_equal(
    get_captures(
      "(RegExr) was (created)",
      "RegExr was created by gskinner.com, and is proudly hosted by Media Temple."
    ),
    c("RegExr was created", "RegExr", "created")
  )
  expect_equal(
    get_captures(
      "RegExr was created",
      "RegExr was created by gskinner.com, and is proudly hosted by Media Temple."
    ),
    c("RegExr was created")
  )
})


test_that("the Design Parser parses design strings correctly", {
  expect_equal(
    design_parser("blue,blue,red"),
    c("blue", "blue", "red")
  )
  expect_equal(
    design_parser("(blue):2,    red"),
    c("blue", "blue", "red")
  )
  expect_equal(
    design_parser("[blue]:2, red"),
    c("blue", "blue", "red")
  )
  expect_equal(
    design_parser("blue,blue,(red):1"),
    c("blue", "blue", "red")
  )
  expect_equal(
    design_parser("blue,blue,[red]:1"),
    c("blue", "blue", "red")
  )
  expect_equal(
    design_parser("[blue,red]:2"),
    c("blue", "blue", "red", "red")
  )
  expect_equal(
    design_parser("(blue,red):2"),
    c("blue", "red", "blue", "red")
  )
  expect_equal(
    design_parser("(blue*,red*):1"),
    c("blue1", "red1")
  )
  expect_equal(
    design_parser("[blue*,red*]:1"),
    c("blue1", "red1")
  )
  expect_equal(
    design_parser("(blue*,red*):1,red3"),
    c("blue4", "red4", "red3")
  )
  expect_equal(
    design_parser("[blue*,red*]:1,red3"),
    c("blue4", "red4", "red3")
  )
  expect_equal(
    design_parser("[blue*,red*]:2"),
    c("blue1", "blue2", "red1", "red2")
  )
  expect_equal(
    design_parser("(blue*,red*):2"),
    c("blue1", "red1", "blue2", "red2")
  )
  expect_equal(
    design_parser("(blue*,red*):2, [blue*,red*]:2"),
    c("blue1", "red1", "blue2", "red2", "blue3", "blue4", "red3", "red4")
  )
})

test_that("The expansion functions ignore the asterisk if told", {
  expect_equal(
    expand_ordered("banana*,papaya*", 2, 3, TRUE),
    c("banana*,banana*,papaya*,papaya*")
  )
  expect_equal(
    expand_intercalated("banana*,papaya*", 2, 3, TRUE),
    c("banana*,papaya*,banana*,papaya*")
  )
})

test_that("split_design() gives expected output and errors", {
  expect_equal(
    split_design(c("blue1", "red1", "blue2", "red2")),
    list(
      groups = c("blue", "red", "blue", "red"),
      pairings = c(1, 1, 2, 2)
    )
  )

  expect_equal(
    split_design(c("blue", "red", "blue", "red")),
    list(
      groups = c("blue", "red", "blue", "red"),
      pairings = c(NA, NA, NA, NA)
    )
  )

  expect_error(
    split_design(c(""))
  )
})
