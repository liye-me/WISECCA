#' Animal-level cross-validation for sparse CCA
#'
#' Selects sparse CCA penalties using leave-one-animal-out
#' cross-validation. Scaling parameters are estimated from training
#' observations only.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#'   A `wise_modules_result` object is also accepted.
#' @param z Numeric matrix or data frame with the same observations as `x`.
#' @param animal Animal identifier for each observation.
#' @param penalty_x_grid Candidate penalties for `x`.
#' @param penalty_z_grid Candidate penalties for `z`.
#' @param max_iter Maximum number of sparse CCA iterations.
#' @param zero_tol Threshold used to identify selected features.
#' @param trace Whether to print fitting progress.
#'
#' @return An object of class `wise_cv_scca_result`.
#'
#' @export
wise_cv_scca <- function(
    x,
    z,
    animal,
    penalty_x_grid = c(0.25, 0.40, 0.55, 0.70, 0.85, 1.00),
    penalty_z_grid = c(0.25, 0.40, 0.55, 0.70, 0.85, 1.00),
    max_iter = 100L,
    zero_tol = 1e-8,
    trace = FALSE
) {
    if (inherits(x, "wise_modules_result")) {
        x <- x$module_eigengenes
    }

    x <- as.matrix(x)
    z <- as.matrix(z)

    if (!is.numeric(x) || !is.numeric(z)) {
        stop("x and z must be numeric matrices.", call. = FALSE)
    }

    if (nrow(x) != nrow(z)) {
        stop("x and z must have the same number of rows.", call. = FALSE)
    }

    if (length(animal) != nrow(x)) {
        stop(
            "animal must have one value for each observation.",
            call. = FALSE
        )
    }

    animal <- as.character(animal)

    if (anyNA(animal) || any(animal == "")) {
        stop("animal contains missing or empty values.", call. = FALSE)
    }

    if (length(unique(animal)) < 3L) {
        stop("At least three animals are required.", call. = FALSE)
    }

    if (
        !is.null(rownames(x)) &&
        !is.null(rownames(z)) &&
        !identical(rownames(x), rownames(z))
    ) {
        stop("The row names of x and z must be identical.", call. = FALSE)
    }

    observation_names <- rownames(x)

    if (is.null(observation_names)) {
        observation_names <- rownames(z)
    }

    if (is.null(observation_names)) {
        observation_names <- paste0("sample_", seq_len(nrow(x)))
    }

    rownames(x) <- observation_names
    rownames(z) <- observation_names

    prepare_grid <- function(grid, feature_count, name) {
        grid <- as.numeric(grid)

        if (
            length(grid) == 0L ||
            any(!is.finite(grid)) ||
            any(grid <= 0) ||
            any(grid > 1)
        ) {
            stop(
                name,
                " must contain values greater than 0 and no greater than 1.",
                call. = FALSE
            )
        }

        lower_bound <- 1 / sqrt(feature_count) + 1e-6

        sort(unique(pmin(1, pmax(grid, lower_bound))))
    }

    penalty_x_grid <- prepare_grid(
        penalty_x_grid,
        ncol(x),
        "penalty_x_grid"
    )

    penalty_z_grid <- prepare_grid(
        penalty_z_grid,
        ncol(z),
        "penalty_z_grid"
    )

    candidates <- expand.grid(
        penalty_x = penalty_x_grid,
        penalty_z = penalty_z_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    folds <- unique(animal)
    oof_store <- vector("list", nrow(candidates))

    cv_results <- data.frame(
        penalty_x = candidates$penalty_x,
        penalty_z = candidates$penalty_z,
        cv_correlation = NA_real_,
        mean_selected_x = NA_real_,
        mean_selected_z = NA_real_,
        stringsAsFactors = FALSE
    )

    for (candidate_index in seq_len(nrow(candidates))) {
        oof_x <- rep(NA_real_, nrow(x))
        oof_z <- rep(NA_real_, nrow(z))

        selected_x <- numeric(length(folds))
        selected_z <- numeric(length(folds))

        reference_weight <- NULL
        successful <- TRUE

        for (fold_index in seq_along(folds)) {
            test_index <- animal == folds[[fold_index]]
            train_index <- !test_index

            fit <- try(
                wise_scca(
                    x = x[train_index, , drop = FALSE],
                    z = z[train_index, , drop = FALSE],
                    penalty_x = candidates$penalty_x[[candidate_index]],
                    penalty_z = candidates$penalty_z[[candidate_index]],
                    n_components = 1L,
                    max_iter = max_iter,
                    standardize = TRUE,
                    trace = trace,
                    zero_tol = zero_tol
                ),
                silent = TRUE
            )

            if (inherits(fit, "try-error")) {
                successful <- FALSE
                break
            }

            weight_x <- fit$weights_x[, 1L]
            weight_z <- fit$weights_z[, 1L]

            x_test <- sweep(
                x[test_index, , drop = FALSE],
                2L,
                fit$center_x,
                "-"
            )

            x_test <- sweep(
                x_test,
                2L,
                fit$scale_x,
                "/"
            )

            z_test <- sweep(
                z[test_index, , drop = FALSE],
                2L,
                fit$center_z,
                "-"
            )

            z_test <- sweep(
                z_test,
                2L,
                fit$scale_z,
                "/"
            )

            score_x <- drop(x_test %*% weight_x)
            score_z <- drop(z_test %*% weight_z)

            if (is.null(reference_weight)) {
                reference_weight <- weight_x
            } else if (sum(reference_weight * weight_x) < 0) {
                score_x <- -score_x
                score_z <- -score_z
                weight_x <- -weight_x
                weight_z <- -weight_z
            }

            oof_x[test_index] <- score_x
            oof_z[test_index] <- score_z

            selected_x[[fold_index]] <- sum(abs(weight_x) > zero_tol)
            selected_z[[fold_index]] <- sum(abs(weight_z) > zero_tol)
        }

        if (
            successful &&
            all(is.finite(oof_x)) &&
            all(is.finite(oof_z))
        ) {
            cv_results$cv_correlation[[candidate_index]] <-
                stats::cor(oof_x, oof_z)

            cv_results$mean_selected_x[[candidate_index]] <-
                mean(selected_x)

            cv_results$mean_selected_z[[candidate_index]] <-
                mean(selected_z)

            oof_store[[candidate_index]] <- list(
                score_x = oof_x,
                score_z = oof_z
            )
        }
    }

    valid <- which(is.finite(cv_results$cv_correlation))

    if (length(valid) == 0L) {
        stop(
            "All cross-validation models failed.",
            call. = FALSE
        )
    }

    selected_total <-
        cv_results$mean_selected_x +
        cv_results$mean_selected_z

    ranking <- order(
        -cv_results$cv_correlation,
        selected_total,
        na.last = TRUE
    )

    best_index <- ranking[[1L]]

    final_model <- wise_scca(
        x = x,
        z = z,
        penalty_x = cv_results$penalty_x[[best_index]],
        penalty_z = cv_results$penalty_z[[best_index]],
        n_components = 1L,
        max_iter = max_iter,
        standardize = TRUE,
        trace = trace,
        zero_tol = zero_tol
    )

    oof_scores <- data.frame(
        observation = observation_names,
        animal = animal,
        score_x = oof_store[[best_index]]$score_x,
        score_z = oof_store[[best_index]]$score_z,
        stringsAsFactors = FALSE
    )

    result <- list(
        best_penalty_x = cv_results$penalty_x[[best_index]],
        best_penalty_z = cv_results$penalty_z[[best_index]],
        cv_correlation = cv_results$cv_correlation[[best_index]],
        cv_results = cv_results,
        oof_scores = oof_scores,
        final_model = final_model,
        folds = folds
    )

    class(result) <- c(
        "wise_cv_scca_result",
        "list"
    )

    result
}
