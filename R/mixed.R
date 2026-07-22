#' Fit mixed-effects models
#'
#' Fits additive and age-interaction mixed-effects models for an
#' imaging outcome and a molecular module. Animal is included as a
#' random intercept.
#'
#' @param outcome Numeric imaging outcome.
#' @param module Numeric molecular module score.
#' @param meta Data frame containing sample metadata.
#' @param animal_col Animal identifier column in `meta`.
#' @param brain_col Brain-region column in `meta`.
#' @param age_col Age or stage column in `meta`.
#' @param standardize Whether to standardize `outcome` and `module`.
#' @param reml Whether to use restricted maximum likelihood.
#' @param reference_age Optional reference level for age.
#' @param adjust Adjustment method for pairwise slope comparisons.
#' @param max_iter Maximum number of optimizer evaluations.
#'
#' @return An object of class `wise_mixed_result`.
#'
#' @export
wise_mixed <- function(
    outcome,
    module,
    meta,
    animal_col = "animal",
    brain_col = "brain",
    age_col = "age",
    standardize = TRUE,
    reml = TRUE,
    reference_age = NULL,
    adjust = "tukey",
    max_iter = 100000L
) {
    if (!is.data.frame(meta)) {
        stop("meta must be a data frame.", call. = FALSE)
    }

    required_cols <- c(animal_col, brain_col, age_col)
    missing_cols <- setdiff(required_cols, colnames(meta))

    if (length(missing_cols) > 0L) {
        stop(
            "meta is missing required columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    n <- nrow(meta)

    if (
        !is.numeric(outcome) ||
        length(outcome) != n ||
        any(!is.finite(outcome))
    ) {
        stop(
            "outcome must be a finite numeric vector with one value per row.",
            call. = FALSE
        )
    }

    if (
        !is.numeric(module) ||
        length(module) != n ||
        any(!is.finite(module))
    ) {
        stop(
            "module must be a finite numeric vector with one value per row.",
            call. = FALSE
        )
    }

    for (column in required_cols) {
        values <- as.character(meta[[column]])

        if (anyNA(values) || any(values == "")) {
            stop(
                column,
                " contains missing or empty values.",
                call. = FALSE
            )
        }
    }

    if (stats::sd(outcome) == 0) {
        stop("outcome must have non-zero variance.", call. = FALSE)
    }

    if (stats::sd(module) == 0) {
        stop("module must have non-zero variance.", call. = FALSE)
    }

    make_factor <- function(x) {
        if (is.factor(x)) {
            return(droplevels(x))
        }

        factor(x, levels = unique(x))
    }

    model_data <- data.frame(
        outcome = as.numeric(outcome),
        module = as.numeric(module),
        age = make_factor(meta[[age_col]]),
        brain = make_factor(meta[[brain_col]]),
        animal = make_factor(meta[[animal_col]]),
        stringsAsFactors = FALSE
    )

    if (nlevels(model_data$animal) < 3L) {
        stop("At least three animals are required.", call. = FALSE)
    }

    if (nlevels(model_data$age) < 2L) {
        stop("At least two age levels are required.", call. = FALSE)
    }

    if (nlevels(model_data$brain) < 2L) {
        stop("At least two brain-region levels are required.", call. = FALSE)
    }

    if (!is.null(reference_age)) {
        if (!reference_age %in% levels(model_data$age)) {
            stop(
                "reference_age is not present in the age variable.",
                call. = FALSE
            )
        }

        model_data$age <- stats::relevel(
            model_data$age,
            ref = reference_age
        )
    }

    outcome_center <- if (standardize) {
        mean(model_data$outcome)
    } else {
        0
    }

    outcome_scale <- if (standardize) {
        stats::sd(model_data$outcome)
    } else {
        1
    }

    module_center <- if (standardize) {
        mean(model_data$module)
    } else {
        0
    }

    module_scale <- if (standardize) {
        stats::sd(model_data$module)
    } else {
        1
    }

    model_data$outcome_z <- (
        model_data$outcome - outcome_center
    ) / outcome_scale

    model_data$module_z <- (
        model_data$module - module_center
    ) / module_scale

    control <- lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = as.integer(max_iter))
    )

    additive_model <- lme4::lmer(
        outcome_z ~ module_z + age + brain + (1 | animal),
        data = model_data,
        REML = reml,
        control = control,
        na.action = stats::na.fail
    )

    interaction_model <- lme4::lmer(
        outcome_z ~ module_z * age + brain + (1 | animal),
        data = model_data,
        REML = reml,
        control = control,
        na.action = stats::na.fail
    )

    extract_fixed <- function(model) {
        coefficient_matrix <- stats::coef(summary(model))

        data.frame(
            term = rownames(coefficient_matrix),
            estimate = coefficient_matrix[, "Estimate"],
            std_error = coefficient_matrix[, "Std. Error"],
            statistic = coefficient_matrix[, "t value"],
            row.names = NULL,
            stringsAsFactors = FALSE
        )
    }

    stage_trends <- emmeans::emtrends(
        interaction_model,
        specs = "age",
        var = "module_z",
        lmer.df = "asymptotic"
    )

    stage_slopes <- as.data.frame(
        summary(
            stage_trends,
            infer = c(TRUE, TRUE)
        )
    )

    slope_contrasts <- as.data.frame(
        summary(
            emmeans::contrast(
                stage_trends,
                method = "pairwise",
                adjust = adjust
            ),
            infer = c(TRUE, TRUE)
        )
    )

    convergence_messages <- function(model) {
        messages <- model@optinfo$conv$lme4$messages

        if (is.null(messages)) {
            return(character(0))
        }

        as.character(messages)
    }

    result <- list(
        additive_model = additive_model,
        interaction_model = interaction_model,
        additive_fixed_effects = extract_fixed(additive_model),
        interaction_fixed_effects = extract_fixed(interaction_model),
        stage_slopes = stage_slopes,
        slope_contrasts = slope_contrasts,
        model_data = model_data,
        diagnostics = list(
            additive_singular = lme4::isSingular(
                additive_model,
                tol = 1e-4
            ),
            interaction_singular = lme4::isSingular(
                interaction_model,
                tol = 1e-4
            ),
            additive_convergence = convergence_messages(
                additive_model
            ),
            interaction_convergence = convergence_messages(
                interaction_model
            )
        ),
        scaling = list(
            outcome_center = outcome_center,
            outcome_scale = outcome_scale,
            module_center = module_center,
            module_scale = module_scale
        ),
        parameters = list(
            standardize = standardize,
            reml = reml,
            reference_age = reference_age,
            adjust = adjust,
            max_iter = as.integer(max_iter)
        )
    )

    class(result) <- c(
        "wise_mixed_result",
        "list"
    )

    result
}
