.wise_call <- function(fun, fixed_args, extra_args) {
    if (!is.list(extra_args)) {
        stop("Additional arguments must be supplied as a list.", call. = FALSE)
    }

    duplicated_args <- intersect(
        names(fixed_args),
        names(extra_args)
    )

    if (length(duplicated_args) > 0L) {
        stop(
            "Reserved arguments cannot be supplied again: ",
            paste(duplicated_args, collapse = ", "),
            call. = FALSE
        )
    }

    do.call(
        fun,
        c(fixed_args, extra_args)
    )
}


#' Run the WISECCA workflow
#'
#' Runs input validation, weighted co-expression module construction,
#' animal-level sparse CCA cross-validation, final sparse CCA fitting,
#' and optional bootstrap, LASSO, and mixed-effects analyses.
#'
#' @param omics Numeric molecular-feature matrix with observations in rows.
#' @param imaging Numeric imaging-feature matrix with observations in rows.
#' @param meta Data frame containing sample metadata.
#' @param id_col Sample identifier column in `meta`.
#' @param animal_col Animal identifier column in `meta`.
#' @param brain_col Brain-region column in `meta`.
#' @param age_col Age or stage column in `meta`.
#' @param power Soft-thresholding power used by WGCNA.
#' @param drop_zero_var Whether zero-variance input features should be removed.
#' @param module_args Additional arguments passed to `wise_modules`.
#' @param cv_args Additional arguments passed to `wise_cv_scca`.
#' @param run_bootstrap Whether to run animal-block bootstrap analysis.
#' @param bootstrap_args Additional arguments passed to `wise_bootstrap`.
#' @param run_lasso Whether to run LASSO analysis.
#' @param lasso_outcomes Imaging features used as LASSO outcomes. If `NULL`,
#'   all imaging features are analyzed.
#' @param lasso_args Additional arguments passed to `wise_lasso`.
#' @param run_mixed Whether to run one mixed-effects model.
#' @param mixed_outcome Imaging feature used as the mixed-model outcome.
#' @param mixed_module Module eigengene used as the mixed-model predictor.
#' @param mixed_args Additional arguments passed to `wise_mixed`.
#'
#' @return An object of class `wisecca_result`.
#'
#' @export
wisecca <- function(
    omics,
    imaging,
    meta,
    id_col = "sample_id",
    animal_col = "animal",
    brain_col = "brain",
    age_col = "age",
    power = 6,
    drop_zero_var = TRUE,
    module_args = list(),
    cv_args = list(),
    run_bootstrap = TRUE,
    bootstrap_args = list(),
    run_lasso = FALSE,
    lasso_outcomes = NULL,
    lasso_args = list(),
    run_mixed = FALSE,
    mixed_outcome = NULL,
    mixed_module = NULL,
    mixed_args = list()
) {
    input <- wise_check_input(
        omics = omics,
        imaging = imaging,
        meta = meta,
        id_col = id_col,
        animal_col = animal_col,
        brain_col = brain_col,
        age_col = age_col,
        drop_zero_var = drop_zero_var
    )

    modules <- .wise_call(
        wise_modules,
        fixed_args = list(
            omics = input$omics,
            power = power
        ),
        extra_args = module_args
    )

    module_matrix <- modules$module_eigengenes
    imaging_matrix <- input$imaging
    animal <- input$meta[[animal_col]]

    cv_scca <- .wise_call(
        wise_cv_scca,
        fixed_args = list(
            x = module_matrix,
            z = imaging_matrix,
            animal = animal
        ),
        extra_args = cv_args
    )

    bootstrap_result <- NULL

    if (isTRUE(run_bootstrap)) {
        bootstrap_result <- .wise_call(
            wise_bootstrap,
            fixed_args = list(
                x = module_matrix,
                z = imaging_matrix,
                animal = animal,
                cv_result = cv_scca
            ),
            extra_args = bootstrap_args
        )
    }

    lasso_result <- NULL

    if (isTRUE(run_lasso)) {
        if (is.null(lasso_outcomes)) {
            lasso_outcomes <- colnames(imaging_matrix)
        }

        invalid_outcomes <- setdiff(
            lasso_outcomes,
            colnames(imaging_matrix)
        )

        if (length(invalid_outcomes) > 0L) {
            stop(
                "Unknown LASSO outcomes: ",
                paste(invalid_outcomes, collapse = ", "),
                call. = FALSE
            )
        }

        lasso_result <- lapply(
            lasso_outcomes,
            function(outcome_name) {
                .wise_call(
                    wise_lasso,
                    fixed_args = list(
                        x = module_matrix,
                        y = imaging_matrix[, outcome_name],
                        animal = animal
                    ),
                    extra_args = lasso_args
                )
            }
        )

        names(lasso_result) <- lasso_outcomes
    }

    mixed_result <- NULL

    if (isTRUE(run_mixed)) {
        if (
            is.null(mixed_outcome) ||
            length(mixed_outcome) != 1L ||
            !mixed_outcome %in% colnames(imaging_matrix)
        ) {
            stop(
                "mixed_outcome must identify one imaging feature.",
                call. = FALSE
            )
        }

        if (
            is.null(mixed_module) ||
            length(mixed_module) != 1L ||
            !mixed_module %in% colnames(module_matrix)
        ) {
            stop(
                "mixed_module must identify one module eigengene.",
                call. = FALSE
            )
        }

        mixed_result <- .wise_call(
            wise_mixed,
            fixed_args = list(
                outcome = imaging_matrix[, mixed_outcome],
                module = module_matrix[, mixed_module],
                meta = input$meta,
                animal_col = animal_col,
                brain_col = brain_col,
                age_col = age_col
            ),
            extra_args = mixed_args
        )
    }

    result <- list(
        input = input,
        modules = modules,
        cv_scca = cv_scca,
        scca = cv_scca$final_model,
        bootstrap = bootstrap_result,
        lasso = lasso_result,
        mixed = mixed_result,
        parameters = list(
            power = power,
            run_bootstrap = run_bootstrap,
            run_lasso = run_lasso,
            run_mixed = run_mixed
        ),
        call = match.call()
    )

    class(result) <- c(
        "wisecca_result",
        "list"
    )

    result
}
