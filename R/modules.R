.wise_with_wgcna_cor <- function(expr) {
    target_env <- globalenv()
    symbols <- c("cor", "bicor")

    existed <- vapply(
        symbols,
        exists,
        logical(1),
        envir = target_env,
        inherits = FALSE
    )

    previous <- stats::setNames(vector("list", length(symbols)), symbols)

    for (symbol in symbols[existed]) {
        previous[[symbol]] <- get(
            symbol,
            envir = target_env,
            inherits = FALSE
        )
    }

    on.exit(
        {
            for (symbol in symbols) {
                if (existed[[symbol]]) {
                    assign(
                        symbol,
                        previous[[symbol]],
                        envir = target_env
                    )
                } else if (
                    exists(
                        symbol,
                        envir = target_env,
                        inherits = FALSE
                    )
                ) {
                    rm(list = symbol, envir = target_env)
                }
            }
        },
        add = TRUE
    )

    assign("cor", WGCNA::cor, envir = target_env)
    assign("bicor", WGCNA::bicor, envir = target_env)

    force(expr)
}


#' Construct weighted co-expression modules
#'
#' Constructs weighted co-expression modules from molecular data using
#' WGCNA and calculates module eigengenes and module membership values.
#'
#' @param omics Numeric matrix or data frame with observations in rows
#'   and molecular features in columns. A `wise_input` object is also accepted.
#' @param power Positive soft-thresholding power.
#' @param network_type Network type. One of `"signed"`, `"unsigned"`,
#'   or `"signed hybrid"`.
#' @param cor_type Correlation type. One of `"pearson"` or `"bicor"`.
#' @param min_module_size Minimum number of features in a module.
#' @param merge_cut_height Height used to merge similar modules.
#' @param deep_split Module splitting sensitivity from 0 to 4.
#' @param max_block_size Maximum number of features in one block.
#' @param random_seed Random seed used by WGCNA.
#' @param n_threads Number of parallel threads. Zero uses automatic detection.
#' @param verbose WGCNA verbosity level.
#'
#' @return An object of class `wise_modules_result` containing module
#'   assignments, module eigengenes, module membership values, removed
#'   zero-variance features, parameters, and the original WGCNA result.
#'
#' @export
wise_modules <- function(
    omics,
    power = 6,
    network_type = c("signed", "unsigned", "signed hybrid"),
    cor_type = c("pearson", "bicor"),
    min_module_size = 30L,
    merge_cut_height = 0.25,
    deep_split = 2L,
    max_block_size = 5000L,
    random_seed = 54321L,
    n_threads = 0L,
    verbose = 0L
) {
    network_type <- match.arg(network_type)
    cor_type <- match.arg(cor_type)

    if (inherits(omics, "wise_input")) {
        omics <- omics$omics
    }

    if (!is.matrix(omics) && !is.data.frame(omics)) {
        stop(
            "omics must be a numeric matrix, data frame, or wise_input object.",
            call. = FALSE
        )
    }

    omics <- as.matrix(omics)

    if (!is.numeric(omics)) {
        stop(
            "omics must contain only numeric values.",
            call. = FALSE
        )
    }

    storage.mode(omics) <- "double"

    if (nrow(omics) < 4L) {
        stop(
            "omics must contain at least four observations.",
            call. = FALSE
        )
    }

    if (ncol(omics) < 4L) {
        stop(
            "omics must contain at least four features.",
            call. = FALSE
        )
    }

    if (
        is.null(colnames(omics)) ||
        anyNA(colnames(omics)) ||
        any(colnames(omics) == "")
    ) {
        stop(
            "omics must have non-empty feature names.",
            call. = FALSE
        )
    }

    if (anyDuplicated(colnames(omics))) {
        stop(
            "omics contains duplicated feature names.",
            call. = FALSE
        )
    }

    if (any(!is.finite(omics))) {
        stop(
            "omics contains NA, NaN, or infinite values.",
            call. = FALSE
        )
    }

    if (is.null(rownames(omics))) {
        rownames(omics) <- paste0(
            "sample_",
            seq_len(nrow(omics))
        )
    }

    if (anyDuplicated(rownames(omics))) {
        stop(
            "omics contains duplicated observation names.",
            call. = FALSE
        )
    }

    if (
        length(power) != 1L ||
        !is.numeric(power) ||
        !is.finite(power) ||
        power <= 0
    ) {
        stop(
            "power must be a positive numeric value.",
            call. = FALSE
        )
    }

    if (
        length(min_module_size) != 1L ||
        !is.numeric(min_module_size) ||
        !is.finite(min_module_size) ||
        min_module_size < 2L
    ) {
        stop(
            "min_module_size must be at least 2.",
            call. = FALSE
        )
    }

    min_module_size <- as.integer(min_module_size)

    if (min_module_size > ncol(omics)) {
        stop(
            "min_module_size cannot exceed the number of features.",
            call. = FALSE
        )
    }

    if (
        length(merge_cut_height) != 1L ||
        !is.numeric(merge_cut_height) ||
        !is.finite(merge_cut_height) ||
        merge_cut_height < 0 ||
        merge_cut_height > 1
    ) {
        stop(
            "merge_cut_height must be between 0 and 1.",
            call. = FALSE
        )
    }

    if (
        length(deep_split) != 1L ||
        !is.numeric(deep_split) ||
        !is.finite(deep_split) ||
        !deep_split %in% 0:4
    ) {
        stop(
            "deep_split must be an integer from 0 to 4.",
            call. = FALSE
        )
    }

    deep_split <- as.integer(deep_split)

    if (
        length(max_block_size) != 1L ||
        !is.numeric(max_block_size) ||
        !is.finite(max_block_size) ||
        max_block_size < min_module_size
    ) {
        stop(
            "max_block_size must be at least min_module_size.",
            call. = FALSE
        )
    }

    max_block_size <- as.integer(max_block_size)

    if (
        length(n_threads) != 1L ||
        !is.numeric(n_threads) ||
        !is.finite(n_threads) ||
        n_threads < 0
    ) {
        stop(
            "n_threads must be a non-negative integer.",
            call. = FALSE
        )
    }

    n_threads <- as.integer(n_threads)

    zero_variance <- vapply(
        seq_len(ncol(omics)),
        function(j) {
            length(unique(omics[, j])) <= 1L
        },
        logical(1)
    )

    removed_features <- colnames(omics)[zero_variance]

    if (any(zero_variance)) {
        warning(
            "Zero-variance features were removed: ",
            paste(removed_features, collapse = ", "),
            call. = FALSE
        )

        omics <- omics[, !zero_variance, drop = FALSE]
    }

    if (ncol(omics) < min_module_size) {
        stop(
            "Too few features remain after removing zero-variance features.",
            call. = FALSE
        )
    }

    tom_type <- if (network_type == "unsigned") {
        "unsigned"
    } else {
        "signed"
    }

    fit <- .wise_with_wgcna_cor(
        WGCNA::blockwiseModules(
            datExpr = omics,
            power = power,
            networkType = network_type,
            corType = cor_type,
            TOMType = tom_type,
            minModuleSize = min_module_size,
            mergeCutHeight = merge_cut_height,
            deepSplit = deep_split,
            maxBlockSize = max_block_size,
            randomSeed = as.integer(random_seed),
            nThreads = n_threads,
            numericLabels = FALSE,
            saveTOMs = FALSE,
            useCorOptionsThroughout = TRUE,
            verbose = as.integer(verbose)
        )
    )

    if (!isTRUE(fit$MEsOK)) {
        stop(
            "WGCNA failed to calculate valid module eigengenes.",
            call. = FALSE
        )
    }

    module_colors <- as.character(fit$colors)
    names(module_colors) <- colnames(omics)

    module_eigengenes <- as.matrix(fit$MEs)
    rownames(module_eigengenes) <- rownames(omics)

    if ("MEgrey" %in% colnames(module_eigengenes)) {
        module_eigengenes <- module_eigengenes[
            ,
            colnames(module_eigengenes) != "MEgrey",
            drop = FALSE
        ]
    }

    if (ncol(module_eigengenes) == 0L) {
        stop(
            "No non-grey modules were detected.",
            call. = FALSE
        )
    }

    module_membership <- if (cor_type == "bicor") {
        WGCNA::bicor(
            omics,
            module_eigengenes,
            use = "pairwise.complete.obs"
        )
    } else {
        WGCNA::cor(
            omics,
            module_eigengenes,
            use = "pairwise.complete.obs"
        )
    }

    module_membership <- as.matrix(module_membership)

    colnames(module_membership) <- sub(
        "^ME",
        "kME",
        colnames(module_membership)
    )

    rownames(module_membership) <- colnames(omics)

    assigned_kme <- vapply(
        seq_along(module_colors),
        function(i) {
            membership_name <- paste0(
                "kME",
                module_colors[[i]]
            )

            if (
                module_colors[[i]] == "grey" ||
                !membership_name %in% colnames(module_membership)
            ) {
                return(NA_real_)
            }

            module_membership[i, membership_name]
        },
        numeric(1)
    )

    module_assignment <- data.frame(
        feature = colnames(omics),
        module = unname(module_colors),
        kME = assigned_kme,
        stringsAsFactors = FALSE
    )

    result <- list(
        module_assignment = module_assignment,
        module_eigengenes = module_eigengenes,
        module_membership = module_membership,
        removed_zero_variance = removed_features,
        parameters = list(
            power = power,
            network_type = network_type,
            cor_type = cor_type,
            min_module_size = min_module_size,
            merge_cut_height = merge_cut_height,
            deep_split = deep_split,
            max_block_size = max_block_size,
            random_seed = as.integer(random_seed),
            n_threads = n_threads
        ),
        wgcna_result = fit
    )

    class(result) <- c(
        "wise_modules_result",
        "list"
    )

    result
}
