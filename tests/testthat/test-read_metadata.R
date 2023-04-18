testthat::test_that("import_msorganizer_xlm: Template with all information", {
    d <- midar::import_msorganizer_xlm("20_MSTemplate_Creator_forTest.xlsm")
    dd <- readRDS("20_MSTemplate_Creator_forTest.rds")
    expect_identical(d, dd)
})

