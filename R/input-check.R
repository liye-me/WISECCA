#' Check WISECCA input data
#'
#' Validates the formats and sample correspondence of omics, imaging,
#' and metadata inputs.
#'
#' @param omics Numeric matrix or data frame with observations in rows
#'   and molecular features in columns.
#' @param imaging Numeric matrix or data frame with observations in rows
#'   and imaging features in columns.
#' @param meta Data frame containing sample metadata.
#' @param id_col Name of the sample identifier column in `meta`.
#' @param animal_col Name of the animal identifier column in `meta`.
#' @param brain_col Name of the brain-region column in `meta`.
#' @param age_col Name of the age or stage column in `meta`.
#' @param drop_zero_var Logical value indicating whether zero-variance
#'   features should be removed.
#'
#' @return An object of class `wise_input`.
#' @export
wise_check_input <- function(
    omics,
    imaging,
    meta,
    id_col = "sample_id",
    animal_col = "animal",
    brain_col = "brain",
    age_col = "age",
    drop_zero_var = TRUE
) {
    to_numeric_matrix <- function(x, name) {
        if (!is.matrix(x) && !is.data.frame(x)) {
            stop(name, " must be a matrix or data frame.", call. = FALSE)
        }

        x <- as.matrix(x)

        if (!is.numeric(x)) {
            stop(name, " must contain only numeric values.", call. = FALSE)
        }

        if (ncol(x) == 0L) {
            stop(name, " must contain at least one feature.", call. = FALSE)
        }

        if (is.null(colnames(x)) || any(colnames(x) == "")) {
            stop(name, " must have non-empty column names.", call. = FALSE)
        }

        if (anyDuplicated(colnames(x))) {
            stop(name, " contains duplicated column names.", call. = FALSE)
        }

        if (any(!is.finite(x))) {
            stop(name, " contains NA, NaN, or infinite values.", call. = FALSE)
        }

        storage.mode(x) <- "double"
        x
    }

    has_custom_rownames <- function(x) {
        rn <- rownames(x)

        !is.null(rn) &&
            !identical(rn, as.character(seq_len(nrow(x))))
    }

    find_zero_variance <- function(x) {
        vapply(
            seq_len(ncol(x)),
            function(j) length(unique(x[, j])) <= 1L,
            logical(1)
        )
    }

    if (!is.data.frame(meta)) {
        stop("meta must be a data frame.", call. = FALSE)
    }

    required_cols <- c(id_col, animal_col, brain_col, age_col)
    missing_cols <- setdiff(required_cols, colnames(meta))

    if (length(missing_cols) > 0L) {
        stop(
            "meta is missing required columns: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
        )
    }

    omics_has_rownames <- has_custom_rownames(omics)
    imaging_has_rownames <- has_custom_rownames(imaging)

    omics <- to_numeric_matrix(omics, "omics")
    imaging <- to_numeric_matrix(imaging, "imaging")

    if (
        nrow(omics) != nrow(imaging) ||
        nrow(omics) != nrow(meta)
    ) {
        stop(
            "omics, imaging, and meta must have the same number of rows.",
            call. = FALSE
        )
    }

    sample_id <- as.character(meta[[id_col]])

    if (anyNA(sample_id) || any(sample_id == "")) {
        stop(
            id_col,
            " contains missing or empty sample identifiers.",
            call. = FALSE
        )
    }

    if (anyDuplicated(sample_id)) {
        stop(
            id_col,
            " contains duplicated sample identifiers.",
            call. = FALSE
        )
    }

    for (column in c(animal_col, brain_col, age_col)) {
        values <- as.character(meta[[column]])

        if (anyNA(values) || any(values == "")) {
            stop(
                column,
                " contains missing or empty values.",
                call. = FALSE
            )
        }
    }

    if (
        omics_has_rownames &&
        !identical(rownames(omics), sample_id)
    ) {
        stop(
            "omics row names do not match sample identifiers in meta.",
            call. = FALSE
        )
    }

    if (
        imaging_has_rownames &&
        !identical(rownames(imaging), sample_id)
    ) {
        stop(
            "imaging row names do not match sample identifiers in meta.",
            call. = FALSE
        )
    }

    rownames(omics) <- sample_id
    rownames(imaging) <- sample_id
    rownames(meta) <- sample_id

    omics_zero <- find_zero_variance(omics)
    imaging_zero <- find_zero_variance(imaging)

    removed <- list(
        omics = colnames(omics)[omics_zero],
        imaging = colnames(imaging)[imaging_zero]
    )

    if ((any(omics_zero) || any(imaging_zero)) && !drop_zero_var) {
        stop("Zero-variance features were detected.", call. = FALSE)
    }

    if (drop_zero_var) {
        omics <- omics[, !omics_zero, drop = FALSE]
        imaging <- imaging[, !imaging_zero, drop = FALSE]
    }

    if (ncol(omics) == 0L) {
        stop(
            "No omics features remain after input validation.",
            call. = FALSE
        )
    }

    if (ncol(imaging) == 0L) {
        stop(
            "No imaging features remain after input validation.",
            call. = FALSE
        )
    }

    result <- list(
        omics = omics,
        imaging = imaging,
        meta = meta,
        removed_zero_variance = removed
    )

    class(result) <- c("wise_input", "list")
    result
}
