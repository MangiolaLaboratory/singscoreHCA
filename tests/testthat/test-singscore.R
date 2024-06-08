
test_that("singscore", {
   singscoreHCA(
     c("CD3G", "CD8A"),
     tissue = "lung",
     cell_type = "t_nk"
    ) |>
    expect_no_error()
})
