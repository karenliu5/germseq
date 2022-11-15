test_that("graph input is correct", {
  EXPECTED_OUTPUT <- c(39, 76, 45, 70, 55, 60)
  graph_output <- get_graph_input1(atlas1006_output)
  graph_val <- graph_output$value

  expect_equal(sum(graph_val == EXPECTED_OUTPUT), 6)
  expect_equal(dim(graph_output), c(6,3))
})


