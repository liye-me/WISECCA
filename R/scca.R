#' Fit sparse canonical correlation analysis
#'
#' Fits sparse canonical correlation analysis between two data matrices
#' using penalized matrix decomposition.
#'
#' @param x Numeric matrix or data frame with observations in rows and
#'   features in columns. A `wise_modules_result` object is also accepted.
#' @param z Numeric matrix or data frame with the same observations as `x`.
#' @param penalty_x Sparsity parameter for `x`. Larger values produce
#'   less sparse canonical weights.
#' @param penalty_z Sparsity parameter for `z`. Larger values produce
#'   less sparse canonical weights.
#' @param n_components Number of canonical components.
#' @param max_iter Maximum number of iterations.
#' @param standardize Whether to center and scale each feature.
#' @param trace Whether to print fitting progress.
#' @param zero_tol Threshold used to identify selected features.
#'
#' @return An object of class `wise_scca_result`.
#'
#' @export
wise_scca <- function(
    x,
    z,
    penalty_x = 0.8,
    penalty_z = 0.8,
    n_components = 1L,
    max_iter = 100L,
    standardize = TRUE,
    trace = FALSE,
    zero_tol = 1e-8
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

        if (any(!is.finite(data))) {
            stop(
                name,
                " contains NA, NaN, or infinite values.",
                call. = FALSE
            )
        }

        storage.mode(data) <- "double"
        data
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

    if (ncol(x) < 2L || ncol(z) < 2L) {
        stop(
            "x and z must each contain at least two features.",
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

    if (is.null(rownames(x)) && is.null(rownames(z))) {
        observation_names <- paste0(
            "sample_",
            seq_len(nrow(x))
        )
    } else if (!is.null(rownames(x))) {
        observation_names <- rownames(x)
    } else {
        observation_names <- rownames(z)
    }

    if (anyDuplicated(observation_names)) {
        stop(
            "Observation names must be unique.",
            call. = FALSE
        )
    }

    rownames(x) <- observation_names
    rownames(z) <- observation_names

    validate_penalty <- function(penalty, data, name) {
        lower_bound <- 1 / sqrt(ncol(data))

        if (
            length(penalty) != 1L ||
            !is.numeric(penalty) ||
            !is.finite(penalty) ||
            penalty < lower_bound ||
            penalty > 1
        ) {
            stop(
                name,
                " must be between ",
                signif(lower_bound, 4),
                " and 1.",
                call. = FALSE
            )
        }

        penalty
    }

    penalty_x <- validate_penalty(
        penalty_x,
        x,
        "penalty_x"
    )

    penalty_z <- validate_penalty(
        penalty_z,
        z,
        "penalty_z"
    )

    max_components <- min(
        ncol(x),
        ncol(z),
        nrow(x) - 1L
    )

    if (
        length(n_components) != 1L ||
        !is.numeric(n_components) ||
        !is.finite(n_components) ||
        n_components != as.integer(n_components) ||
        n_components < 1L ||
        n_components > max_components
    ) {
        stop(
            "n_components must be an integer between 1 and ",
            max_components,
            ".",
            call. = FALSE
        )
    }

    n_components <- as.integer(n_components)

    if (
        length(max_iter) != 1L ||
        !is.numeric(max_iter) ||
        !is.finite(max_iter) ||
        max_iter != as.integer(max_iter) ||
        max_iter < 1L
    ) {
        stop(
            "max_iter must be a positive integer.",
            call. = FALSE
        )
    }

    max_iter <- as.integer(max_iter)

    scale_data <- function(data, standardize) {
        if (standardize) {
            center <- colMeans(data)
            scale_value <- apply(data, 2L, stats::sd)

            if (
                any(!is.finite(scale_value)) ||
                any(scale_value <= 0)
            ) {
                stop(
                    "Zero-variance features were detected.",
                    call. = FALSE
                )
            }

            scaled <- sweep(data, 2L, center, "-")
            scaled <- sweep(scaled, 2L, scale_value, "/")
        } else {
            center <- rep(0, ncol(data))
            scale_value <- rep(1, ncol(data))
            scaled <- data
        }

        names(center) <- colnames(data)
        names(scale_value) <- colnames(data)

        list(
            data = scaled,
            center = center,
            scale = scale_value
        )
    }

    x_info <- scale_data(x, standardize)
    z_info <- scale_data(z, standardize)

    fit <- PMA::CCA(
        x = x_info$data,
        z = z_info$data,
        typex = "standard",
        typez = "standard",
        penaltyx = penalty_x,
        penaltyz = penalty_z,
        K = n_components,
        niter = max_iter,
        trace = trace,
        standardize = FALSE
    )

    weights_x <- as.matrix(fit$u)
    weights_z <- as.matrix(fit$v)

    component_names <- paste0(
        "CC",
        seq_len(n_components)
    )

    rownames(weights_x) <- colnames(x)
    rownames(weights_z) <- colnames(z)
    colnames(weights_x) <- component_names
    colnames(weights_z) <- component_names

    scores_x <- x_info$data %*% weights_x
    scores_z <- z_info$data %*% weights_z

    rownames(scores_x) <- observation_names
    rownames(scores_z) <- observation_names
    colnames(scores_x) <- component_names
    colnames(scores_z) <- component_names

    for (component in seq_len(n_components)) {
        reference_index <- which.max(
            abs(weights_x[, component])
        )

        if (weights_x[reference_index, component] < 0) {
            weights_x[, component] <- -weights_x[, component]
            weights_z[, component] <- -weights_z[, component]
            scores_x[, component] <- -scores_x[, component]
            scores_z[, component] <- -scores_z[, component]
        }
    }

    canonical_correlations <- vapply(
        seq_len(n_components),
        function(component) {
            stats::cor(
                scores_x[, component],
                scores_z[, component]
            )
        },
        numeric(1)
    )

    names(canonical_correlations) <- component_names

    selected_x <- lapply(
        seq_len(n_components),
        function(component) {
            rownames(weights_x)[
                abs(weights_x[, component]) > zero_tol
            ]
        }
    )

    selected_z <- lapply(
        seq_len(n_components),
        function(component) {
            rownames(weights_z)[
                abs(weights_z[, component]) > zero_tol
            ]
        }
    )

    names(selected_x) <- component_names
    names(selected_z) <- component_names

    result <- list(
        weights_x = weights_x,
        weights_z = weights_z,
        scores_x = scores_x,
        scores_z = scores_z,
        canonical_correlations = canonical_correlations,
        selected_x = selected_x,
        selected_z = selected_z,
        center_x = x_info$center,
        scale_x = x_info$scale,
        center_z = z_info$center,
        scale_z = z_info$scale,
        parameters = list(
            penalty_x = penalty_x,
            penalty_z = penalty_z,
            n_components = n_components,
            max_iter = max_iter,
            standardize = standardize,
            zero_tol = zero_tol
        ),
        pma_result = fit
    )

    class(result) <- c(
        "wise_scca_result",
        "list"
    )

    result
}
