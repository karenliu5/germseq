test_that("graph input is correct", {
  EXPECTED_OUTPUT_val <- c(52, 18, 14, 31)
  EXPECTED_OUTPUT_names <- c("no methods", "one method",
                             "two methods", "all methods")

  graph_output <- get_graph_input2(atlas1006_output)
  graph_val <- graph_output$value
  graph_names <- graph_output$method

  expect_equal(sum(graph_val == EXPECTED_OUTPUT_val), 4)
  expect_equal(sum(graph_names == EXPECTED_OUTPUT_names), 4)
  expect_equal(dim(graph_output), c(4,2))
})



