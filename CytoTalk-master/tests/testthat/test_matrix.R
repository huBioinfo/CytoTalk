test_that("grouping is identical to extraction", {
    c <- 6

    mat <- matrix(runif(10 * c), ncol = c)
    meta <- paste(as.integer((seq_len(c) - 1) / 2))

    lst <- CytoTalk:::new_named_list(mat, meta)

    expect_identical(
        CytoTalk::extract_group(meta[1], lst),
        CytoTalk::group_meta(lst)[[meta[1]]]
    )
})
