#' Animal-block bootstrap for sparse CCA
#'
#' Assesses sparse CCA stability by resampling animals with replacement.
#' All observations from the same animal are resampled as one block.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#'   A `wise_modules_result` object is also accepted.
#' @param z Numeric matrix or data frame with the same observations as `x`.
#' @param animal Animal identifier for each observation.
#' @param penalty_x Sparsity parameter for `x`.
#' @param penalty_z Sparsity parameter for `z`.
#' @param cv_result Optional `wise_cv_scca_result` object used to obtain
#'   the selected penalties.
#' @param n_boot Number of bootstrap replicates.
#' @param selection_threshold Minimum selection frequency for stability.
#' @param sign_threshold Minimum sign consistency for stability.
#' @param max_iter Maximum number of sparse CCA iterations.
#' @param zero_tol Threshold used to identify non-zero weights.
#' @param seed Random seed.
#' @param trace Whether to print progress.
#'
#' @return An object of class `wise_bootstrap_result`.
#'
#' @export
wise_bootstrap <- function(
    x,
    z,
    animal,
    penalty_x = NULL,
    penalty_z = NULL,
    cv_result = NULL,
    n_boot = 500L,
    selection_threshold = 0.70,
    sign_threshold = 0.80,
    max_iter = 100L,
    zero_tol = 1e-8,
    seed = 123L,
    trace = FALSE
) {
    if (inherits(x, "wise_modules_result")) {
        x <- x$module_eigengenes
    }

    to_matrix <- function(data, name) {
        if (!is.matrix(data) && !is.data.frame(data)) {
            stop(
                name,
                " must be a numeric matrix or data frame.",
                call. = FALSE
            )
        }

        data <- as.matrix(data)

        if (!is.numeric(data)) {
            stop(
                name,
                " must contain only numeric values.",
                call. = FALSE
            )
        }

        if (any(!is.finite(data))) {
            stop(
                name,
                " contains NA, NaN, or infinite values.",
                call. = FALSE
            )
        }

        if (
            is.null(colnames(data)) ||
            anyNA(colnames(data)) ||
            any(colnames(data) == "")
        ) {
            stop(
                name,
                " must have non-empty feature names.",
                call. = FALSE
            )
        }

        if (anyDuplicated(colnames(data))) {
            stop(
                name,
                " contains duplicated feature names.",
                call. = FALSE
            )
        }

        storage.mode(data) <- "double"
        data
    }

    validate_probability <- function(value, name) {
        if (
            length(value) != 1L ||
            !is.numeric(value) ||
            !is.finite(value) ||
            value < 0 ||
            value > 1
        ) {
            stop(
                name,
                " must be between 0 and 1.",
                call. = FALSE
            )
        }

        value
    }

    x <- to_matrix(x, "x")
    z <- to_matrix(z, "z")

    if (nrow(x) != nrow(z)) {
        stop(
            "x and z must have the same number of rows.",
            call. = FALSE
        )
    }

    if (nrow(x) < 3L) {
        stop(
            "At least three observations are required.",
            call. = FALSE
        )
    }

    if (
        !is.null(rownames(x)) &&
        !is.null(rownames(z)) &&
        !identical(rownames(x), rownames(z))
    ) {
        stop(
            "The row names of x and z must be identical.",
            call. = FALSE
        )
    }

    if (
        length(animal) != nrow(x) ||
        anyNA(animal) ||
        any(as.character(animal) == "")
    ) {
        stop(
            "animal must contain one non-missing value per observation.",
            call. = FALSE
        )
    }

    animal <- as.character(animal)
    animals <- unique(animal)

    if (length(animals) < 3L) {
        stop(
            "At least three animals are required.",
            call. = FALSE
        )
    }

    if (!is.null(cv_result)) {
        if (!inherits(cv_result, "wise_cv_scca_result")) {
            stop(
                "cv_result must be a wise_cv_scca_result object.",
                call. = FALSE
            )
        }

        penalty_x <- cv_result$best_penalty_x
        penalty_z <- cv_result$best_penalty_z
    }

    if (is.null(penalty_x) || is.null(penalty_z)) {
        stop(
            "penalty_x and penalty_z must be supplied when cv_result is NULL.",
            call. = FALSE
        )
    }

    if (
        length(n_boot) != 1L ||
        !is.numeric(n_boot) ||
        !is.finite(n_boot) ||
        n_boot != as.integer(n_boot) ||
        n_boot < 10L
    ) {
        stop(
            "n_boot must be an integer of at least 10.",
            call. = FALSE
        )
    }

    n_boot <- as.integer(n_boot)

    selection_threshold <- validate_probability(
        selection_threshold,
        "selection_threshold"
    )

    sign_threshold <- validate_probability(
        sign_threshold,
        "sign_threshold"
    )

    observation_names <- rownames(x)

    if (is.null(observation_names)) {
        observation_names <- rownames(z)
    }

    if (is.null(observation_names)) {
        observation_names <- paste0(
            "sample_",
            seq_len(nrow(x))
        )
    }

    rownames(x) <- observation_names
    rownames(z) <- observation_names

    reference_model <- wise_scca(
        x = x,
        z = z,
        penalty_x = penalty_x,
        penalty_z = penalty_z,
        n_components = 1L,
        max_iter = max_iter,
        standardize = TRUE,
        trace = FALSE,
        zero_tol = zero_tol
    )

    reference_x <- reference_model$weights_x[, 1L]
    reference_z <- reference_model$weights_z[, 1L]

    weights_x <- matrix(
        NA_real_,
        nrow = n_boot,
        ncol = ncol(x),
        dimnames = list(
            paste0("bootstrap_", seq_len(n_boot)),
            colnames(x)
        )
    )

    weights_z <- matrix(
        NA_real_,
        nrow = n_boot,
        ncol = ncol(z),
        dimnames = list(
            paste0("bootstrap_", seq_len(n_boot)),
            colnames(z)
        )
    )

    canonical_correlations <- rep(NA_real_, n_boot)
    errors <- rep(NA_character_, n_boot)

    set.seed(as.integer(seed))

    for (bootstrap_index in seq_len(n_boot)) {
        sampled_animals <- sample(
            animals,
            size = length(animals),
            replace = TRUE
        )

        sampled_rows <- unlist(
            lapply(
                sampled_animals,
                function(current_animal) {
                    which(animal == current_animal)
                }
            ),
            use.names = FALSE
        )

        x_boot <- x[sampled_rows, , drop = FALSE]
        z_boot <- z[sampled_rows, , drop = FALSE]

        bootstrap_names <- paste0(
            "bootstrap_",
            bootstrap_index,
            "_row_",
            seq_len(nrow(x_boot))
        )

        rownames(x_boot) <- bootstrap_names
        rownames(z_boot) <- bootstrap_names

        fit <- try(
            wise_scca(
                x = x_boot,
                z = z_boot,
                penalty_x = penalty_x,
                penalty_z = penalty_z,
                n_components = 1L,
                max_iter = max_iter,
                standardize = TRUE,
                trace = FALSE,
                zero_tol = zero_tol
            ),
            silent = TRUE
        )

        if (inherits(fit, "try-error")) {
            errors[[bootstrap_index]] <- as.character(fit)
            next
        }

        current_x <- fit$weights_x[, 1L]
        current_z <- fit$weights_z[, 1L]

        alignment <- sum(
            reference_x * current_x,
            na.rm = TRUE
        )

        if (abs(alignment) <= zero_tol) {
            alignment <- sum(
                reference_z * current_z,
                na.rm = TRUE
            )
        }

        if (is.finite(alignment) && alignment < 0) {
            current_x <- -current_x
            current_z <- -current_z
        }

        weights_x[bootstrap_index, ] <- current_x
        weights_z[bootstrap_index, ] <- current_z

        canonical_correlations[[bootstrap_index]] <-
            fit$canonical_correlations[[1L]]

        if (isTRUE(trace)) {
            message(
                "Completed bootstrap replicate ",
                bootstrap_index,
                " of ",
                n_boot,
                "."
            )
        }
    }

    successful <- is.finite(canonical_correlations)
    n_successful <- sum(successful)

    if (n_successful == 0L) {
        stop(
            "All bootstrap models failed.",
            call. = FALSE
        )
    }

    summarize_weights <- function(
        weight_matrix,
        reference_weight
    ) {
        weight_matrix <- weight_matrix[
            successful,
            ,
            drop = FALSE
        ]

        selected <- abs(weight_matrix) > zero_tol

        selection_frequency <- colMeans(selected)

        positive_count <- colSums(
            weight_matrix > zero_tol
        )

        negative_count <- colSums(
            weight_matrix < -zero_tol
        )

        selected_count <- positive_count + negative_count

        sign_consistency <- ifelse(
            selected_count > 0L,
            pmax(positive_count, negative_count) / selected_count,
            NA_real_
        )

        dominant_sign <- ifelse(
            selected_count == 0L,
            NA_character_,
            ifelse(
                positive_count >= negative_count,
                "positive",
                "negative"
            )
        )

        mean_weight <- colMeans(weight_matrix)

        median_weight <- apply(
            weight_matrix,
            2L,
            stats::median
        )

        ci_lower <- apply(
            weight_matrix,
            2L,
            stats::quantile,
            probs = 0.025,
            names = FALSE
        )

        ci_upper <- apply(
            weight_matrix,
            2L,
            stats::quantile,
            probs = 0.975,
            names = FALSE
        )

        stable <- (
            selection_frequency >= selection_threshold &
                sign_consistency >= sign_threshold
        )

        data.frame(
            feature = colnames(weight_matrix),
            reference_weight = unname(reference_weight),
            mean_weight = unname(mean_weight),
            median_weight = unname(median_weight),
            ci_lower = unname(ci_lower),
            ci_upper = unname(ci_upper),
            selection_frequency = unname(selection_frequency),
            sign_consistency = unname(sign_consistency),
            dominant_sign = unname(dominant_sign),
            stable = unname(stable),
            stringsAsFactors = FALSE
        )
    }

    stability_x <- summarize_weights(
        weights_x,
        reference_x
    )

    stability_z <- summarize_weights(
        weights_z,
        reference_z
    )

    bootstrap_correlations <- data.frame(
        replicate = seq_len(n_boot),
        correlation = canonical_correlations,
        successful = successful,
        error = errors,
        stringsAsFactors = FALSE
    )

    result <- list(
        reference_model = reference_model,
        stability_x = stability_x,
        stability_z = stability_z,
        stable_x = stability_x$feature[stability_x$stable],
        stable_z = stability_z$feature[stability_z$stable],
        weights_x = weights_x,
        weights_z = weights_z,
        bootstrap_correlations = bootstrap_correlations,
        n_successful = n_successful,
        n_failed = n_boot - n_successful,
        parameters = list(
            penalty_x = penalty_x,
            penalty_z = penalty_z,
            n_boot = n_boot,
            selection_threshold = selection_threshold,
            sign_threshold = sign_threshold,
            max_iter = max_iter,
            zero_tol = zero_tol,
            seed = as.integer(seed)
        )
    )

    class(result) <- c(
        "wise_bootstrap_result",
        "list"
    )

    result
}
