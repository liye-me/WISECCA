#' Fit animal-level cross-validated LASSO
#'
#' Fits a Gaussian LASSO model using animal-level cross-validation.
#' All observations from the same animal are assigned to the same fold.
#'
#' @param x Numeric predictor matrix with observations in rows and
#'   features in columns. A `wise_modules_result` object is also accepted.
#' @param y Numeric outcome vector.
#' @param animal Animal identifier for each observation.
#' @param lambda_choice Lambda selection rule. Either `"lambda.1se"`
#'   or `"lambda.min"`.
#' @param standardize Whether predictors should be standardized.
#' @param intercept Whether the model should include an intercept.
#' @param max_iter Maximum number of coordinate-descent iterations.
#' @param zero_tol Threshold used to identify selected features.
#'
#' @return An object of class `wise_lasso_result`.
#'
#' @export
wise_lasso <- function(
    x,
    y,
    animal,
    lambda_choice = c("lambda.1se", "lambda.min"),
    standardize = TRUE,
    intercept = TRUE,
    max_iter = 100000L,
    zero_tol = 1e-8
) {
    lambda_choice <- match.arg(lambda_choice)

    if (inherits(x, "wise_modules_result")) {
        x <- x$module_eigengenes
    }

    if (!is.matrix(x) && !is.data.frame(x)) {
        stop(
            "x must be a numeric matrix or data frame.",
            call. = FALSE
        )
    }

    x <- as.matrix(x)

    if (!is.numeric(x)) {
        stop(
            "x must contain only numeric values.",
            call. = FALSE
        )
    }

    if (any(!is.finite(x))) {
        stop(
            "x contains NA, NaN, or infinite values.",
            call. = FALSE
        )
    }

    if (
        is.null(colnames(x)) ||
        anyNA(colnames(x)) ||
        any(colnames(x) == "")
    ) {
        stop(
            "x must have non-empty feature names.",
            call. = FALSE
        )
    }

    if (anyDuplicated(colnames(x))) {
        stop(
            "x contains duplicated feature names.",
            call. = FALSE
        )
    }

    storage.mode(x) <- "double"

    if (
        !is.numeric(y) ||
        length(y) != nrow(x) ||
        any(!is.finite(y))
    ) {
        stop(
            "y must be a finite numeric vector with one value per observation.",
            call. = FALSE
        )
    }

    y <- as.numeric(y)

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

    if (stats::sd(y) == 0) {
        stop(
            "y must have non-zero variance.",
            call. = FALSE
        )
    }

    zero_variance <- vapply(
        seq_len(ncol(x)),
        function(j) stats::sd(x[, j]) == 0,
        logical(1)
    )

    removed_features <- colnames(x)[zero_variance]

    if (any(zero_variance)) {
        warning(
            "Zero-variance features were removed: ",
            paste(removed_features, collapse = ", "),
            call. = FALSE
        )

        x <- x[, !zero_variance, drop = FALSE]
    }

    if (ncol(x) == 0L) {
        stop(
            "No predictors remain after removing zero-variance features.",
            call. = FALSE
        )
    }

    observation_names <- rownames(x)

    if (is.null(observation_names)) {
        observation_names <- paste0(
            "sample_",
            seq_len(nrow(x))
        )
    }

    if (anyDuplicated(observation_names)) {
        stop(
            "Observation names must be unique.",
            call. = FALSE
        )
    }

    rownames(x) <- observation_names

    fold_id <- match(animal, animals)

    cv_fit <- glmnet::cv.glmnet(
        x = x,
        y = y,
        family = "gaussian",
        alpha = 1,
        foldid = fold_id,
        nfolds = length(animals),
        type.measure = "mse",
        standardize = standardize,
        intercept = intercept,
        keep = TRUE,
        maxit = as.integer(max_iter)
    )

    selected_lambda <- cv_fit[[lambda_choice]]

    coefficient_matrix <- as.matrix(
        stats::coef(
            cv_fit,
            s = lambda_choice
        )
    )

    coefficient_vector <- coefficient_matrix[, 1L]
    names(coefficient_vector) <- rownames(coefficient_matrix)

    intercept_value <- unname(
        coefficient_vector["(Intercept)"]
    )

    feature_coefficients <- coefficient_vector[
        names(coefficient_vector) != "(Intercept)"
    ]

    selected_features <- names(feature_coefficients)[
        abs(feature_coefficients) > zero_tol
    ]

    lambda_index <- which.min(
        abs(cv_fit$lambda - selected_lambda)
    )

    oof_prediction <- drop(
        cv_fit$fit.preval[, lambda_index]
    )

    valid_prediction <- is.finite(oof_prediction)

    if (sum(valid_prediction) < 3L) {
        stop(
            "Cross-validated predictions could not be calculated.",
            call. = FALSE
        )
    }

    cv_correlation <- stats::cor(
        y[valid_prediction],
        oof_prediction[valid_prediction]
    )

    cv_rmse <- sqrt(
        mean(
            (
                y[valid_prediction] -
                    oof_prediction[valid_prediction]
            )^2
        )
    )

    coefficients <- data.frame(
        feature = names(feature_coefficients),
        coefficient = unname(feature_coefficients),
        selected = abs(feature_coefficients) > zero_tol,
        stringsAsFactors = FALSE
    )

    predictions <- data.frame(
        observation = observation_names,
        animal = animal,
        observed = y,
        predicted = oof_prediction,
        residual = y - oof_prediction,
        stringsAsFactors = FALSE
    )

    result <- list(
        selected_lambda = selected_lambda,
        lambda_choice = lambda_choice,
        intercept = intercept_value,
        coefficients = coefficients,
        selected_features = selected_features,
        predictions = predictions,
        cv_correlation = cv_correlation,
        cv_rmse = cv_rmse,
        fold_id = fold_id,
        removed_zero_variance = removed_features,
        cv_model = cv_fit
    )

    class(result) <- c(
        "wise_lasso_result",
        "list"
    )

    result
}
