library(microbiome)
data(atlas1006)

# Tests for error checking
test_that("threshold check works", {

  # compare_DAA_methods()
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = -1),
               "prevThr")
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 1),
               "prevThr")
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 0),
               "two values")
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 0.5),
               "two values")
})


test_that("String check works", {

  # compare_DAA_methods()
  expect_error(compare_DAA_methods(atlas1006, 100, prevThr = 0.5),
               "string")
  expect_error(compare_DAA_methods(atlas1006, atlas1006, prevThr = 0.5),
               "string")
  expect_error(compare_DAA_methods(atlas1006, c("a", "b"), prevThr = 0.5),
               "string")
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 0.5),
               "two values")
})


test_that("Valid group check works", {

  # compare_DAA_methods()
  expect_error(compare_DAA_methods(atlas1006, "asdf", prevThr = 0.5),
               "found")
  expect_error(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 0.5),
               "two values")
})


test_that("Check that it successfully removes NAs", {

  atlas1006_narm <- subset_samples(atlas1006, ! is.na(bmi_group))

  # compare_DAA_methods()
  expect_warning(compare_DAA_methods(atlas1006, "bmi_group", prevThr = 0.5),
               "missing values")
  expect_error(compare_DAA_methods(atlas1006_narm, "bmi_group", prevThr = 0.5),
               "two_values")

})


# Test for content

test_that("Check prevalence filter", {
  atlas1006_clean <- phyloseq::subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
  prevdf <- microbiome::prevalence(atlas1006_clean) # Get prevalence

  prevThr1 <- 0.1
  mask_1 <- as.logical(prevdf > prevThr1)
  expected_size_1 <- sum(mask_1)

  prevThr2 <- 0.5
  mask_2 <- as.logical(prevdf > prevThr2)
  expected_size_2 <- sum(mask_2)

  output1 <- compare_DAA_methods(atlas1006_clean, "bmi_group", prevThr = prevThr1)
  output2 <- compare_DAA_methods(atlas1006_clean, "bmi_group", prevThr = prevThr2)

  # compare_DAA_methods()
  expect_equal(nrow(output1), expected_size_1)
  expect_equal(nrow(output2), expected_size_2)
})


