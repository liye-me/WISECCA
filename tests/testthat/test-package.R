test_that("public functions are exported", {
    expected_exports <- c(
        "wise_check_input",
        "wise_modules",
        "wise_scca",
        "wise_cv_scca",
        "wise_bootstrap",
        "wise_lasso",
        "wise_mixed",
        "wisecca"
    )

    package_exports <- getNamespaceExports("WISECCA")

    expect_true(
        all(expected_exports %in% package_exports)
    )
})

test_that("plot method is registered", {
    plot_method <- getS3method(
        "plot",
        "wisecca_result",
        optional = TRUE
    )

    expect_true(is.function(plot_method))
})
