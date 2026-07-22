#' Plot WISECCA results
#'
#' Creates diagnostic and result plots for a `wisecca_result` object.
#'
#' @param x A `wisecca_result` object.
#' @param type Plot type. One of `"scores"`, `"cv"`, `"weights"`,
#'   `"bootstrap"`, `"lasso"`, or `"mixed"`.
#' @param component Canonical component to plot.
#' @param side Data side to plot. Either `"x"` or `"z"`.
#' @param top_n Maximum number of features to display.
#' @param outcome LASSO outcome name. The first available outcome is used
#'   when this argument is `NULL`.
#' @param ... Additional arguments reserved for future use.
#'
#' @return A `ggplot` object.
#'
#' @export
plot.wisecca_result <- function(
    x,
    type = c(
        "scores",
        "cv",
        "weights",
        "bootstrap",
        "lasso",
        "mixed"
    ),
    component = 1L,
    side = c("x", "z"),
    top_n = 20L,
    outcome = NULL,
    ...
) {
    if (!inherits(x, "wisecca_result")) {
        stop(
            "x must be a wisecca_result object.",
            call. = FALSE
        )
    }

    type <- match.arg(type)
    side <- match.arg(side)

    if (
        length(component) != 1L ||
        !is.numeric(component) ||
        component < 1L ||
        component != as.integer(component)
    ) {
        stop(
            "component must be a positive integer.",
            call. = FALSE
        )
    }

    component <- as.integer(component)

    if (
        length(top_n) != 1L ||
        !is.numeric(top_n) ||
        top_n < 1L ||
        top_n != as.integer(top_n)
    ) {
        stop(
            "top_n must be a positive integer.",
            call. = FALSE
        )
    }

    top_n <- as.integer(top_n)

    if (type == "scores") {
        scores_x <- x$scca$scores_x
        scores_z <- x$scca$scores_z

        if (
            component > ncol(scores_x) ||
            component > ncol(scores_z)
        ) {
            stop(
                "The requested component is not available.",
                call. = FALSE
            )
        }

        correlation <- x$scca$canonical_correlations[[component]]

        plot_data <- data.frame(
            score_x = scores_x[, component],
            score_z = scores_z[, component]
        )

        return(
            ggplot2::ggplot(
                plot_data,
                ggplot2::aes(x = score_x, y = score_z)
            ) +
                ggplot2::geom_point(size = 2) +
                ggplot2::geom_smooth(
                    method = "lm",
                    formula = y ~ x,
                    se = FALSE
                ) +
                ggplot2::labs(
                    x = paste0("X score, CC", component),
                    y = paste0("Z score, CC", component),
                    title = "Sparse CCA canonical scores",
                    subtitle = paste0(
                        "Correlation = ",
                        signif(correlation, 3)
                    )
                ) +
                ggplot2::theme_bw()
        )
    }

    if (type == "cv") {
        plot_data <- x$cv_scca$cv_results
        plot_data <- plot_data[
            is.finite(plot_data$cv_correlation),
            ,
            drop = FALSE
        ]

        if (nrow(plot_data) == 0L) {
            stop(
                "No valid cross-validation results are available.",
                call. = FALSE
            )
        }

        plot_data$penalty_x <- factor(
            plot_data$penalty_x,
            levels = sort(unique(plot_data$penalty_x))
        )

        plot_data$penalty_z <- factor(
            plot_data$penalty_z,
            levels = sort(unique(plot_data$penalty_z))
        )

        plot_data$label <- sprintf(
            "%.2f",
            plot_data$cv_correlation
        )

        return(
            ggplot2::ggplot(
                plot_data,
                ggplot2::aes(
                    x = penalty_x,
                    y = penalty_z,
                    fill = cv_correlation
                )
            ) +
                ggplot2::geom_tile() +
                ggplot2::geom_text(
                    ggplot2::aes(label = label)
                ) +
                ggplot2::labs(
                    x = "Penalty X",
                    y = "Penalty Z",
                    fill = "CV correlation",
                    title = "Sparse CCA penalty selection"
                ) +
                ggplot2::theme_bw()
        )
    }

    if (type == "weights") {
        weights <- if (side == "x") {
            x$scca$weights_x
        } else {
            x$scca$weights_z
        }

        if (component > ncol(weights)) {
            stop(
                "The requested component is not available.",
                call. = FALSE
            )
        }

        plot_data <- data.frame(
            feature = rownames(weights),
            weight = weights[, component],
            stringsAsFactors = FALSE
        )

        plot_data <- plot_data[
            order(abs(plot_data$weight), decreasing = TRUE),
            ,
            drop = FALSE
        ]

        plot_data <- utils::head(plot_data, top_n)

        plot_data$feature <- factor(
            plot_data$feature,
            levels = rev(plot_data$feature)
        )

        return(
            ggplot2::ggplot(
                plot_data,
                ggplot2::aes(x = feature, y = weight)
            ) +
                ggplot2::geom_col() +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    x = NULL,
                    y = "Canonical weight",
                    title = paste0(
                        "Sparse CCA weights: ",
                        toupper(side),
                        " side"
                    )
                ) +
                ggplot2::theme_bw()
        )
    }

    if (type == "bootstrap") {
        if (is.null(x$bootstrap)) {
            stop(
                "Bootstrap results are not available.",
                call. = FALSE
            )
        }

        plot_data <- if (side == "x") {
            x$bootstrap$stability_x
        } else {
            x$bootstrap$stability_z
        }

        plot_data <- plot_data[
            order(
                plot_data$selection_frequency,
                decreasing = TRUE
            ),
            ,
            drop = FALSE
        ]

        plot_data <- utils::head(plot_data, top_n)

        plot_data$feature <- factor(
            plot_data$feature,
            levels = rev(plot_data$feature)
        )

        threshold <- x$bootstrap$parameters$selection_threshold

        return(
            ggplot2::ggplot(
                plot_data,
                ggplot2::aes(
                    x = feature,
                    y = selection_frequency,
                    fill = stable
                )
            ) +
                ggplot2::geom_col() +
                ggplot2::geom_hline(
                    yintercept = threshold,
                    linetype = 2
                ) +
                ggplot2::coord_flip() +
                ggplot2::labs(
                    x = NULL,
                    y = "Selection frequency",
                    fill = "Stable",
                    title = paste0(
                        "Bootstrap stability: ",
                        toupper(side),
                        " side"
                    )
                ) +
                ggplot2::theme_bw()
        )
    }

    if (type == "lasso") {
        if (is.null(x$lasso) || length(x$lasso) == 0L) {
            stop(
                "LASSO results are not available.",
                call. = FALSE
            )
        }

        if (is.null(outcome)) {
            outcome <- names(x$lasso)[[1L]]
        }

        if (!outcome %in% names(x$lasso)) {
            stop(
                "The requested LASSO outcome is not available.",
                call. = FALSE
            )
        }

        plot_data <- x$lasso[[outcome]]$predictions

        return(
            ggplot2::ggplot(
                plot_data,
                ggplot2::aes(
                    x = observed,
                    y = predicted
                )
            ) +
                ggplot2::geom_point(size = 2) +
                ggplot2::geom_abline(
                    slope = 1,
                    intercept = 0,
                    linetype = 2
                ) +
                ggplot2::labs(
                    x = "Observed",
                    y = "Cross-validated prediction",
                    title = paste0(
                        "LASSO prediction: ",
                        outcome
                    ),
                    subtitle = paste0(
                        "CV correlation = ",
                        signif(
                            x$lasso[[outcome]]$cv_correlation,
                            3
                        )
                    )
                ) +
                ggplot2::theme_bw()
        )
    }

    if (type == "mixed") {
        if (is.null(x$mixed)) {
            stop(
                "Mixed-effects results are not available.",
                call. = FALSE
            )
        }

        slopes <- x$mixed$stage_slopes

        trend_column <- grep(
            "\\.trend$",
            colnames(slopes),
            value = TRUE
        )

        if (length(trend_column) == 0L) {
            stop(
                "The mixed-model trend column was not found.",
                call. = FALSE
            )
        }

        trend_column <- trend_column[[1L]]

        lower_column <- intersect(
            c("asymp.LCL", "lower.CL"),
            colnames(slopes)
        )

        upper_column <- intersect(
            c("asymp.UCL", "upper.CL"),
            colnames(slopes)
        )

        plot_data <- data.frame(
            age = as.character(slopes$age),
            estimate = slopes[[trend_column]],
            stringsAsFactors = FALSE
        )

        plot_data$age <- factor(
            plot_data$age,
            levels = unique(plot_data$age)
        )

        if (
            length(lower_column) > 0L &&
            length(upper_column) > 0L
        ) {
            plot_data$lower <- slopes[[lower_column[[1L]]]]
            plot_data$upper <- slopes[[upper_column[[1L]]]]
        }

        plot_object <- ggplot2::ggplot(
            plot_data,
            ggplot2::aes(
                x = age,
                y = estimate,
                group = 1
            )
        ) +
            ggplot2::geom_hline(
                yintercept = 0,
                linetype = 2
            ) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 2) +
            ggplot2::labs(
                x = "Age",
                y = "Estimated module slope",
                title = "Age-specific mixed-model slopes"
            ) +
            ggplot2::theme_bw()

        if (
            "lower" %in% colnames(plot_data) &&
            "upper" %in% colnames(plot_data)
        ) {
            plot_object <- plot_object +
                ggplot2::geom_errorbar(
                    ggplot2::aes(
                        ymin = lower,
                        ymax = upper
                    ),
                    width = 0.15
                )
        }

        return(plot_object)
    }

    stop("Unsupported plot type.", call. = FALSE)
}
