#' Run paired DESeq2 and DEXSeq analyses
#' 
#' With \code{paired_diff} you can run a paired differential gene expression and
#' splicing analysis. The function expects a counts matrix or a
#' \code{\link[SummarizedExperiment:SummarizedExperiment]{SummarizedExperiment}}
#' or
#' \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}} object as input.
#' A preliminary prefiltering step is performed to remove genes with a summed 
#' count lower than the provided threshold. Likewise, genes with counts in 
#' only one sample are removed. This step is mostly to speed up differential 
#' analyses, as \code{\link[DESeq2:DESeq]{DESeq2}} will do a stricter filtering.
#' Surrogate Variable Analysis is recommended to allow the analyses to take
#' batch effects etc. into account.
#' After the two differential analyses, the transcript-level p-values will be
#' aggregated to gene-level to allow subsequent Gene-Set Enrichment Analysis.
#' Transcript-level results can be extracted by setting
#' \code{store_results = TRUE}.
#' 
#' 
#' 
#' @inheritParams DESeq2::DESeq
#' @param object A data object of the types matrix, 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}},
#' or \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}}.
#' If a matrix is used, please also provide metadata.
#' @param metadata (Default: \code{NULL}) A metadata file or
#' \code{data.frame} object
#' @param group_col The metadata column specifying the what group each
#' sample is associated with
#' @param sample_col The column in the metadata that specifies the sample IDs
#' (should correspond to column names in \code{object}).
#' Set to \code{"rownames"} if the rownames should be used.
#' @param baseline Group value of baseline samples
#' @param case Group value of case samples
#' @param covariates Name of column(s) in the \code{metadata} that indicate(s)
#' covariates. E.g., c("gender", "tissue_type")
#' @param experiment_title Title of your experiment. Your results will be
#' stored in \code{paste0("results/", experiment_title, "_pairedGSEA.RDS")}.
#' @param run_sva (Default: \code{TRUE})
#' A logical stating whether SVA should be run.
#' @param sva_n (Default: \code{NULL}) Number of surrogate variables to use.
#' \code{NULL} estimates the number from the data. A non-negative integer
#' requests an exact number, subject to the residual degrees-of-freedom guard.
#' @param sva_seed (Default: \code{20260625L}) Seed used by the permutation-based
#' automatic \code{n.sv} estimator.
#' @param sva_max_n (Default: \code{NULL}) Optional maximum for automatic SVA.
#' When supplied, the selected number is also capped at
#' \code{floor(sqrt(n_samples))}. This argument is not used with fixed
#' \code{sva_n}.
#' @param sva_min_residual_df (Default: \code{3L}) Minimum residual degrees of
#' freedom that must remain after adding surrogate variables.
#' @param use_limma (Default: \code{FALSE})
#' A logical determining if \code{limma+voom} or
#' \code{DESeq2} + \code{DEXSeq} should be used for the analysis
#' @param prefilter (Default: \code{10}) The prefilter threshold,
#' where \code{rowSums}
#' lower than the prefilter threshold will be removed from the count matrix.
#' Set to 0 or \code{FALSE} to prevent prefiltering
#' @param fit_type (Default: \code{"local"}) Either 
#' \code{"parametric", "local", "mean", or "glmGamPoi"}
#' for the type of fitting of dispersions to the mean intensity.
#' @param store_results (Default: \code{FALSE})
#' A logical indicating if results should
#' be stored in the folder \code{"results/"}.
#' @param quiet (Default: \code{FALSE}) Whether to print messages
#' @param expression_only (Default: \code{FALSE}) 
#' A logical that indicates whether to only
#' run \code{\link[DESeq2:DESeq]{DESeq2}} analysis. Not generally recommended.
#' The setting was implemented to make the SVA impact analysis easier
#' @param custom_design (Default: \code{FALSE}) A logical or formula.
#' Can be used to apply a custom design formula for the analysis.
#' Generally not recommended, 
#' as \code{pairedGSEA} will make its own design formula
#' from the group and \code{covariate} columns
#' @param parallel (Default: \code{FALSE}) If FALSE, no parallelization.
#' If TRUE, parallel execution using
#' \code{\link[BiocParallel:bpparam]{BiocParallel}}, see next argument
#' \code{BPPARAM}.
#' @param BPPARAM (Default: 
#' \code{\link[BiocParallel:bpparam]{bpparam()}})
#' An optional
#' parameter object passed internally to
#' \code{\link[BiocParallel:bplapply]{bplapply}}
#' when \code{parallel = TRUE}.
#' If not specified, the parameters last registered with register will be used.
#' @param ... Additional parameters passed to
#' \code{\link[DESeq2:DESeq]{DESeq()}}
#' @family paired
#' @importFrom methods is
#' @importFrom S4Vectors DataFrame
#' @import SummarizedExperiment
#' @import DESeq2
#' @import DEXSeq
#' @import BiocParallel
#' @return A DFrame of aggregated pvalues. SVA selection diagnostics are stored
#' in \code{S4Vectors::metadata(result)$pairedGSEA_sva}.
#' @usage 
#' paired_diff(
#'     object,
#'     group_col,
#'     sample_col,
#'     baseline,
#'     case,
#'     metadata = NULL,
#'     covariates = NULL,
#'     experiment_title = NULL,
#'     store_results = FALSE,
#'     run_sva = TRUE,
#'     sva_n = NULL,
#'     sva_seed = 20260625L,
#'     sva_max_n = NULL,
#'     sva_min_residual_df = 3L,
#'     use_limma = FALSE,
#'     prefilter = 10,
#'     test = "LRT",
#'     fit_type = "local",
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     expression_only = FALSE,
#'     custom_design = FALSE,
#'     ...
#'     )
#' @examples 
#' 
#' # Run analysis on included example data
#' data("example_se")
#' 
#' diff_results <- paired_diff(
#'     object = example_se[1:15, ],
#'     group_col = "group_nr",
#'     sample_col = "id",
#'     baseline = 1,
#'     case = 2,
#'     experiment_title = "Example",
#'     store_results = FALSE 
#' )
#' 
#' @export
paired_diff <- function(
        object,
        group_col,
        sample_col,
        baseline,
        case,
        metadata = NULL,
        covariates = NULL,
        experiment_title = NULL,
        store_results = FALSE,
        run_sva = TRUE,
        sva_n = NULL,
        sva_seed = 20260625L,
        sva_max_n = NULL,
        sva_min_residual_df = 3L,
        use_limma = FALSE,
        prefilter = 10,
        test = "LRT",
        fit_type = "local",
        quiet = FALSE,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        expression_only = FALSE,
        custom_design = FALSE,
        ...){

    ## Initial error checks
    stopifnot(is(
        c(quiet, use_limma, expression_only, run_sva, store_results, parallel),
        "logical"))
    stopifnot(is(custom_design, "logical") | is(custom_design, "formula"))
    stopifnot(
        "You cannot use 'expression_only' and 'use_limma" =
            !(use_limma & expression_only)
        )
    
    stopifnot(
        "Please provide a valid data type: matrix, SummarizedExperiment,
        DESeqDataSet" = any(class(object) %in% c(
            "matrix", "SummarizedExperiment","DESeqDataSet")))
    
    stopifnot(
        "Cannot store results if experiment title haven't been given" =
            (store_results & !is.null(experiment_title)) | !store_results)
    
    stopifnot(
        "Please provide metadata with your count matrix" =
            (is(object, "matrix") & !is.null(metadata)) | !is(object, "matrix"))
    
    ## Checking column names
    stopifnot(
        "Covariate names must not contain spaces" =
            !grepl(" ", covariates))
    stopifnot(
        "group_col name must not contain spaces" =
            !grepl(" ", group_col))
    stopifnot(
        "sample_col name must not contain spaces" =
            !grepl(" ", sample_col)
        | sample_col == "rownames")

    if(!quiet) message("Running ", experiment_title)

    ## Define design formula
    if(custom_design != FALSE){
        if(custom_design == TRUE & is(object, "DESeqDataSet")) {
            design <- DESeq2::design(object)
        } else if(is(custom_design, "formula")) {
            stopifnot(
                "Please ensure the group_col is in the custom design" =
                    grepl(group_col, as.character(custom_design)[2]))
            design <- custom_design
            
            if(is(object, "DESeqDataSet")){
                if(DESeq2::design(object) != design){
                    message(
                        "OBS: your design will be overwritten to: ",
                        as.character(design))
                }
            }
        }
    } else{
        design <- formularise_vector(c(group_col, covariates))
    }

    ## Convert se to dds
    if (is(object, "SummarizedExperiment")) {
        object <- DESeq2::DESeqDataSet(object, design)
    }

    ## Load metadata
    if(!quiet) message("Preparing metadata")
    if(!(is(object, "matrix")) & is.null(metadata))
        metadata <- SummarizedExperiment::colData(object)

    if(sample_col == "rownames"){
        metadata$sample <- rownames(metadata)
        sample_col <- "sample"
    }
    
    metadata <- prepare_metadata(metadata, group_col, paste(c(baseline, case)))
    
    # ensure columns are factors
    cols_to_factor <- c(group_col, covariates)
    metadata[cols_to_factor] <- lapply(metadata[cols_to_factor], factor)
    
    ## Subsample metadata to only include samples present in the count matrix
    metadata <- metadata[metadata[[sample_col]] %in% colnames(object), ]
    stopifnot(
        "Please ensure that the sample IDs in the metadata matches the
        column names of the count matrix." = nrow(metadata) > 0)
    
    ## Check for presence of undesired characters
    stopifnot(
    "OBS! Some or all sample names contain a '-' or ' ',\
    which will cause downstream methods to complain.\
    Please rename these." = !any(grepl("[- ]", metadata[[sample_col]])))

    ## Check sample_col is in metadata
    stopifnot(
        "Sample column not in metadata" = sample_col %in% colnames(metadata))

    # Add metadata to DESeqDataSet
    if(is(object, "DESeqDataSet")){
        object <- 
            object[, colnames(object) %in% rownames(metadata), drop = FALSE]
        SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(metadata)
    }

    ## Convert count matrix to DESeqDataSet
    if(is(object, "matrix")){
        ## Ensure rows in metadata matches columns in the count matrix
        object <- object[, metadata[[sample_col]]]
        
        ## Create DDS from count matrix
        if(!quiet) message("Converting count matrix to DESeqDataSet")
        object <- convert_matrix_to_dds(object, metadata, design)
        }

    stopifnot(
        "Please ensure the rownames have the format 'gene:transcript'" =
            (expression_only | (!expression_only & 
                                    grepl(":", rownames(object)[1]))))

    # Rename object variable to dds
    dds <- object; rm(object)

    # Prefiltering
    if(prefilter) dds <- pre_filter(dds, prefilter, quiet = quiet)

    # Detect surrogate variables and retain an audit record for every mode
    dds <- run_sva(
        dds,
        group_col = group_col,
        enabled = run_sva,
        sva_n = sva_n,
        sva_seed = sva_seed,
        sva_max_n = sva_max_n,
        sva_min_residual_df = sva_min_residual_df,
        quiet = quiet)
    
    if(use_limma){
        # Run limma
        results <- run_limma(
            dds = dds, 
            group_col = group_col,
            baseline = baseline,
            case = case,
            experiment_title = experiment_title,
            store_results = store_results,
            quiet = quiet
        )
        # Extract results
        expression_results <- results$expression
        splicing_results <- results$splicing
        rm(results)
    } else{ # DESeq2 + DEXSeq
        # Run DESeq2
        expression_results <- run_deseq(
            dds,
            group_col = group_col,
            baseline = baseline,
            case = case,
            fit_type = fit_type,
            test = test,
            experiment_title = experiment_title,
            store_results = store_results,
            quiet = quiet,
            parallel = parallel,
            BPPARAM = BPPARAM,
            ...
        )
        
        ## If only DGE is requested (used for SVA evaluations in the paper)
        if(expression_only){
            expression_aggregated <- aggregate_pvalue(
                expression_results, type = "expression")
            expression_aggregated <- attach_sva_diagnostics(
                expression_aggregated, dds)
            
            if(store_results) store_result(
                expression_aggregated, paste0(
                    experiment_title,
                    "_aggregated_pvals.RDS"),
                "gene pvalue aggregation")
            return(expression_aggregated)
        }
        
        # Run DEXSeq
        splicing_results <- run_dexseq(
            dds,
            group_col = group_col,
            baseline = baseline,
            case = case,
            experiment_title = experiment_title,
            store_results = store_results,
            quiet = quiet,
            parallel = parallel,
            BPPARAM = BPPARAM
        )
    }

    
    # Aggregate p values
    if(!quiet) {
        if(!is(experiment_title, "NULL"))
            message(experiment_title, " is analysed.")
        message("Aggregating p values")
    }

    expression_aggregated <- aggregate_pvalue(
        expression_results, type = "expression")
    
    splicing_aggregated <- aggregate_pvalue(
        splicing_results, type = "splicing")

    
    aggregated_pvals <- merge(
        expression_aggregated,
        splicing_aggregated,
        by = "gene",
        suffixes = c("_expression", "_splicing"),
        all = TRUE
    )


    aggregated_pvals <- attach_sva_diagnostics(aggregated_pvals, dds)

    if(store_results) store_result(
        aggregated_pvals, paste0(experiment_title, "_aggregated_pvals.RDS"),
        "gene pvalue aggregation")

    if(!quiet) message("Done.")


    return(aggregated_pvals)
}


#' Prepare metadata
#' 
#' This internal function reads a filepath for a metadata file
#' (or a data.frame object), loads it, 
#' and filters out the rows that do not contain data from the
#' relevant patient group, as defined in the \code{comparison} argument.
#' 
#' @examples 
#' \dontrun{
#' metadata <- prepare_metadata(metadata = "path/to/metadata.xlsx",
#'     group_col = "group", comparison = "2v1")
#' }
#' 
#' @inheritParams paired_diff
#' @param baseline_case A character vector with baseline and case values
#' @note Suggested: importFrom readxl read_excel importFrom readr read_csv
#' @keywords internal
#' @noRd
#' @return A data.frame
prepare_metadata <- function(metadata, group_col, baseline_case){

    if (is.character(metadata)) {
        if (endsWith(metadata, ".xlsx")) {
            check_missing_package(package = "readxl")
            metadata <- readxl::read_excel(metadata)
        } else if (endsWith(metadata, ".csv")) {
            check_missing_package(package = "readr")
            metadata <- readr::read_csv(metadata)
        }
    } else if (!is(metadata, "data.frame") & !is(metadata, "DataFrame")){
            stop(
            "Please provide path to a metadata file or a data.frame /
            DataFrame object.")}

    if(!(group_col %in% colnames(metadata)))
        stop("Could not find column ", group_col, " in metadata.")
    
    # Remove irrelevant groups
    metadata <- metadata[
        as.character(metadata[[group_col]]) %in% baseline_case,]
    
    # Check both baseline and case is in metadata
    stopifnot(
        "Check for misspellings in baseline or case IDs." = 
            all(baseline_case %in% metadata[[group_col]])
        )

    # Add comparison levels to metadata
    metadata[[group_col]] <- factor(
        metadata[[group_col]], levels = baseline_case)

    return(metadata)
}


#' Convert matrix to \code{\link[DESeq2]{DESeqDataSet}}
#' 
#' This internal function converts a matrix and metadata object into a
#' \code{\link[DESeq2]{DESeqDataSet}}, using the \code{group_col}
#' as the basic design.
#' 
#' @param tx_count Count matrix to convert to dds
#' @inheritParams paired_diff
#' @keywords internal
#' @noRd
#' @return A \code{\link[DESeq2]{DESeqDataSet}}
convert_matrix_to_dds <- function(tx_count, metadata, design){

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = tx_count,
        colData = metadata,
        design = design)
    return(dds)
}

#' Run SVA on \code{\link[DESeq2]{DESeqDataSet}}
#' 
#' This internal function runs a surrogate variable analysis on the
#' count matrix in the \code{\link[DESeq2]{DESeqDataSet}} 
#' The found surrogate variables will then be added to the metadata
#' and the design formula in the \code{\link[DESeq2]{DESeqDataSet}}
#' object to be used in
#' the DGE and DTU analyses.
#'  
#' @inheritParams paired_diff
#' @param dds A DESeqDataSet. See \code{\link[DESeq2]{DESeqDataSet}} 
#' for more information about the object type.
#' @importFrom stats model.matrix cor na.omit
#' @importFrom sva sva
#' @keywords internal
#' @noRd
#' @return A DESeqDataSet
run_sva <- function(
        dds,
        group_col = NULL,
        enabled = TRUE,
        sva_n = NULL,
        sva_seed = 20260625L,
        sva_max_n = NULL,
        sva_min_residual_df = 3L,
        quiet = FALSE) {
    validate_sva_control(
        enabled, sva_n, sva_seed, sva_max_n, sva_min_residual_df)

    metadata <- SummarizedExperiment::colData(dds)
    design <- DESeq2::design(dds)
    design_terms <- attr(stats::terms(design), "term.labels")
    if(is.null(group_col)) group_col <- design_terms[1]
    reduced <- reduce_formula(design, term = group_col)
    mod1 <- stats::model.matrix(design, data = metadata)
    mod0 <- stats::model.matrix(reduced, data = metadata)
    n_samples <- nrow(mod1)
    model_rank <- qr(mod1)$rank
    residual_df_cap <- max(
        0L, n_samples - model_rank - as.integer(sva_min_residual_df))

    if(!enabled) {
        info <- make_sva_diagnostics(
            mode = "disabled",
            n_samples = n_samples,
            model_rank = model_rank,
            min_residual_df = sva_min_residual_df,
            residual_df_cap = residual_df_cap)
        return(set_sva_diagnostics(dds, info))
    }

    if (!quiet) message("Running SVA")
    normalized_counts <- SummarizedExperiment::assay(
        DESeq2::normTransform(dds))

    estimated_n <- NA_integer_
    if(is.null(sva_n)) {
        estimated_n <- estimate_sva_n(
            normalized_counts, mod1, seed = sva_seed)
    }
    selection <- select_sva_n(
        n_samples = n_samples,
        model_rank = model_rank,
        estimated_n = estimated_n,
        sva_n = sva_n,
        sva_max_n = sva_max_n,
        sva_min_residual_df = sva_min_residual_df)

    if(!quiet) {
        if(is.na(selection$estimated_n)) {
            message("Using ", selection$selected_n,
                    " explicitly requested surrogate variables")
        } else{
            message(
                "Estimated ", selection$estimated_n,
                " surrogate variables; selected ", selection$selected_n)
        }
    }

    if(selection$selected_n == 0) {
        info <- make_sva_diagnostics(
            mode = selection$mode,
            seed = if(is.null(sva_n)) sva_seed else NA_integer_,
            requested_n = sva_n,
            estimated_n = selection$estimated_n,
            max_n = sva_max_n,
            square_root_cap = selection$square_root_cap,
            n_samples = n_samples,
            model_rank = model_rank,
            min_residual_df = sva_min_residual_df,
            residual_df_cap = selection$residual_df_cap)
        return(set_sva_diagnostics(dds, info))
    }

    svseq <- sva::sva(
        normalized_counts,
        mod = mod1,
        mod0 = mod0,
        n.sv = selection$selected_n)

    svs <- data.frame(svseq$sv)
    colnames(svs) <- paste0("sv", seq_len(selection$selected_n))
    filtered <- filter_confounding_svs(svs, mod1)
    svs <- filtered$svs
    if(filtered$removed_n > 0 && !quiet) {
        message(
            "Removing ", filtered$removed_n,
            " surrogate variables confounding with the design model.")
    }

    if(ncol(svs) > 0) {
        SummarizedExperiment::colData(dds) <- cbind(metadata, svs)
        DESeq2::design(dds) <- formularise_vector(
            c(design_terms, colnames(svs)))
    }

    final_model <- stats::model.matrix(
        DESeq2::design(dds),
        data = SummarizedExperiment::colData(dds))
    final_model_rank <- qr(final_model)$rank
    retained_n <- ncol(svs)
    info <- make_sva_diagnostics(
        mode = selection$mode,
        seed = if(is.null(sva_n)) sva_seed else NA_integer_,
        requested_n = sva_n,
        estimated_n = selection$estimated_n,
        max_n = sva_max_n,
        square_root_cap = selection$square_root_cap,
        n_samples = n_samples,
        model_rank = model_rank,
        final_model_rank = final_model_rank,
        min_residual_df = sva_min_residual_df,
        residual_df_cap = selection$residual_df_cap,
        selected_n = selection$selected_n,
        retained_n = retained_n,
        removed_confounding_n = filtered$removed_n)

    if(!quiet) {
        message(
            "Retained ", retained_n, " surrogate variables; ",
            info$final_residual_df, " residual degrees of freedom remain\n")
    }
    set_sva_diagnostics(dds, info)
}


validate_sva_control <- function(
        enabled, sva_n, sva_seed, sva_max_n, sva_min_residual_df) {
    if(!is.logical(enabled) || length(enabled) != 1 || is.na(enabled)) {
        stop("'run_sva' must be TRUE or FALSE.", call. = FALSE)
    }

    check_integer <- function(value, name, minimum = 0L, allow_null = TRUE) {
        if(is.null(value) && allow_null) return(invisible(TRUE))
        valid <- is.numeric(value) && length(value) == 1 &&
            !is.na(value) && is.finite(value) &&
            value == floor(value) && value >= minimum &&
            value <= .Machine$integer.max
        if(!valid) {
            stop(
                "'", name, "' must be ",
                if(allow_null) "NULL or " else "",
                "a single integer >= ", minimum, ".", call. = FALSE)
        }
        invisible(TRUE)
    }

    check_integer(sva_n, "sva_n")
    check_integer(sva_seed, "sva_seed", allow_null = FALSE)
    check_integer(sva_max_n, "sva_max_n")
    check_integer(
        sva_min_residual_df, "sva_min_residual_df",
        minimum = 1L, allow_null = FALSE)

    if(!enabled && (!is.null(sva_n) || !is.null(sva_max_n))) {
        stop(
            "'sva_n' and 'sva_max_n' cannot be set when run_sva = FALSE.",
            call. = FALSE)
    }
    if(!is.null(sva_n) && !is.null(sva_max_n)) {
        stop(
            "'sva_max_n' applies only to automatic SVA; do not combine it ",
            "with a fixed 'sva_n'.", call. = FALSE)
    }
    invisible(TRUE)
}


select_sva_n <- function(
        n_samples,
        model_rank,
        estimated_n = NA_integer_,
        sva_n = NULL,
        sva_max_n = NULL,
        sva_min_residual_df = 3L) {
    residual_df_cap <- max(
        0L,
        as.integer(n_samples) - as.integer(model_rank) -
            as.integer(sva_min_residual_df))

    if(!is.null(sva_n)) {
        requested_n <- as.integer(sva_n)
        if(requested_n > residual_df_cap) {
            stop(
                "Requested sva_n = ", requested_n, " would leave fewer than ",
                sva_min_residual_df, " residual degrees of freedom. ",
                "The largest safe value is ", residual_df_cap, ".",
                call. = FALSE)
        }
        return(list(
            mode = "fixed",
            estimated_n = NA_integer_,
            square_root_cap = NA_integer_,
            residual_df_cap = residual_df_cap,
            selected_n = requested_n))
    }

    if(length(estimated_n) != 1 || is.na(estimated_n) ||
            !is.finite(estimated_n) || estimated_n < 0) {
        stop("Automatic SVA requires a valid estimated n.sv.", call. = FALSE)
    }
    estimated_n <- as.integer(estimated_n)

    if(is.null(sva_max_n)) {
        mode <- "automatic"
        square_root_cap <- NA_integer_
        selected_n <- min(estimated_n, residual_df_cap)
    } else{
        mode <- "automatic_capped"
        square_root_cap <- as.integer(floor(sqrt(n_samples)))
        selected_n <- min(
            estimated_n,
            as.integer(sva_max_n),
            square_root_cap,
            residual_df_cap)
    }

    list(
        mode = mode,
        estimated_n = estimated_n,
        square_root_cap = square_root_cap,
        residual_df_cap = residual_df_cap,
        selected_n = as.integer(selected_n))
}


preserve_rng_state <- function(code) {
    seed_name <- ".Random.seed"
    had_seed <- exists(seed_name, envir = .GlobalEnv, inherits = FALSE)
    if(had_seed) {
        old_seed <- get(seed_name, envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
        if(had_seed) {
            assign(seed_name, old_seed, envir = .GlobalEnv)
        } else if(exists(seed_name, envir = .GlobalEnv, inherits = FALSE)) {
            rm(list = seed_name, envir = .GlobalEnv)
        }
    }, add = TRUE)
    code()
}


estimate_sva_n <- function(dat, mod, seed = 20260625L) {
    as.integer(preserve_rng_state(function() {
        sva::num.sv(dat, mod, method = "be", seed = as.integer(seed))
    }))
}


filter_confounding_svs <- function(svs, mod, threshold = 0.8) {
    svs <- as.data.frame(svs)
    if(ncol(svs) == 0 || ncol(mod) <= 1) {
        return(list(svs = svs, removed_n = 0L))
    }

    cors <- suppressWarnings(stats::cor(
        svs, mod[, seq.int(2, ncol(mod)), drop = FALSE]))
    remove <- apply(cors, 1, function(x) {
        any(!is.finite(x) | abs(x) > threshold)
    })
    list(
        svs = svs[, !remove, drop = FALSE],
        removed_n = as.integer(sum(remove)))
}


make_sva_diagnostics <- function(
        mode,
        seed = NA_integer_,
        requested_n = NULL,
        estimated_n = NA_integer_,
        max_n = NULL,
        square_root_cap = NA_integer_,
        n_samples,
        model_rank,
        final_model_rank = model_rank,
        min_residual_df,
        residual_df_cap,
        selected_n = 0L,
        retained_n = 0L,
        removed_confounding_n = 0L) {
    list(
        mode = mode,
        seed = as.integer(seed),
        requested_n = if(is.null(requested_n)) NA_integer_ else
            as.integer(requested_n),
        estimated_n = as.integer(estimated_n),
        max_n = if(is.null(max_n)) NA_integer_ else as.integer(max_n),
        square_root_cap = as.integer(square_root_cap),
        residual_df_cap = as.integer(residual_df_cap),
        selected_n = as.integer(selected_n),
        capped_n = as.integer(selected_n),
        retained_n = as.integer(retained_n),
        removed_confounding_n = as.integer(removed_confounding_n),
        n_samples = as.integer(n_samples),
        model_rank = as.integer(model_rank),
        final_model_rank = as.integer(final_model_rank),
        min_residual_df = as.integer(min_residual_df),
        final_residual_df = as.integer(n_samples - final_model_rank))
}


set_sva_diagnostics <- function(dds, info) {
    object_metadata <- S4Vectors::metadata(dds)
    object_metadata$pairedGSEA_sva <- info
    S4Vectors::metadata(dds) <- object_metadata
    dds
}


attach_sva_diagnostics <- function(result, dds) {
    info <- S4Vectors::metadata(dds)$pairedGSEA_sva
    result_metadata <- S4Vectors::metadata(result)
    result_metadata$pairedGSEA_sva <- info
    S4Vectors::metadata(result) <- result_metadata
    result
}



normalize_differential_test <- function(test) {
    if(!is.character(test) || length(test) != 1 || is.na(test)) {
        stop("'test' must be either \"LRT\" or \"Wald\".", call. = FALSE)
    }
    normalized <- toupper(test)
    if(!(normalized %in% c("LRT", "WALD"))) {
        stop("'test' must be either \"LRT\" or \"Wald\".", call. = FALSE)
    }
    if(normalized == "WALD") "Wald" else "LRT"
}


#' Run DESeq2 analysis
#' 
#' This internal function runs a differential gene expression analysis using
#' DESeq2 (See \code{\link[DESeq2:DESeq]{DESeq}} for more detailed information).
#' Here, the surrogate variables found by \code{\link{run_sva}},
#' if any, will be added to the \code{\link[DESeq2:DESeqDataSet]{DESeqDataSet}}
#' before running the analysis.
#'   
#' @inheritParams paired_diff
#' @inheritParams run_sva
#' @param ... Additional parameters passed to
#' \code{\link[DESeq2:DESeq]{DESeq()}}
#' @keywords internal
#' @noRd
#' @return A DFrame
#' @usage
#' run_deseq(
#'     dds,
#'     group_col,
#'     baseline,
#'     case,
#'     test = "LRT",
#'     fit_type = "local",
#'     experiment_title = "Experiment Title",
#'     store_results = FALSE,
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     ...
#'     )
#' 
run_deseq <- function(
        dds,
        group_col,
        baseline,
        case,
        test = "LRT",
        fit_type = "local",
        experiment_title = "Experiment Title",
        store_results = FALSE,
        quiet = FALSE,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        ...){
    
    if(!quiet) message("Running DESeq2")
    test <- normalize_differential_test(test)
    if(!(group_col %in% all.vars(DESeq2::design(dds)))) {
        stop(
            group_col,
            " should be the name of a factor in the colData of the ",
            "DESeqDataSet", call. = FALSE)
    }
    
    if(test == "LRT") {
        # Reduced model: remove the group term, keep surrogate variables
        # and covariates
        reduced <- reduce_formula(
            DESeq2::design(dds), term = group_col)
        dds <- DESeq2::DESeq(
            dds,
            reduced = reduced,
            test = test,
            parallel = parallel,
            BPPARAM = BPPARAM,
            fitType = fit_type,
            quiet = quiet,
            ...)
    } else{
        dds <- DESeq2::DESeq(
            dds,
            test = test,
            parallel = parallel,
            BPPARAM = BPPARAM,
            fitType = fit_type,
            quiet = quiet,
            ...)
    }
    
    # Store DESeqDataSet with DESeq2 analysis
    if(store_results) {
        store_result(
            dds, paste0(experiment_title, "_dds.RDS"),
            "DESeqDataSet", quiet = quiet)
    }
    
    if(!quiet) message("Extracting results")
    expression_results <- DESeq2::results(
        dds,
        contrast = c(group_col, paste(c(case, baseline))),
        parallel = parallel,
        BPPARAM = BPPARAM)
    
    # Store result extraction
    if(store_results) {
        store_result(
            expression_results,
            paste0(experiment_title, "_expression_results.RDS"),
            "differential expression results", quiet = quiet)
    }
    
    # Separate gene and transcript
    gene_tx <- do.call("rbind", strsplit(
        rownames(expression_results), ":", fixed = TRUE))
    expression_results$gene <- gene_tx[, 1]
    expression_results$transcript <- gene_tx[, 2]
    expression_results$lfc <- expression_results$log2FoldChange
    
    # Reorder columns
    expression_results <- expression_results[, c(
        "gene", "transcript", "lfc", "pvalue", "padj", "baseMean")]
    
    return(expression_results)
}



#' Run DEXSeq analysis
#' 
#' This internal function runs a differential transcript usage analysis using
#' DEXSeq (See \code{\link[DEXSeq:DEXSeq]{DEXSeq}}
#' for more detailed information).
#' Here, the surrogate variables found by \code{\link{run_sva}}, if any,
#' will be added to the \code{\link[DEXSeq:DEXSeqDataSet]{DEXSeqDataSet}}
#' before running the analysis.
#' 
#' @inheritParams paired_diff
#' @inheritParams run_sva
#' @keywords internal
#' @noRd
#' @return A DFrame
#' @usage
#' run_dexseq(
#'     dds,
#'     group_col,
#'     baseline,
#'     case,
#'     experiment_title = "Experiment Title",
#'     store_results = FALSE,
#'     quiet = FALSE,
#'     parallel = FALSE,
#'     BPPARAM = BiocParallel::bpparam()
#'     )
#'     
run_dexseq <- function(
        dds,
        group_col,
        baseline,
        case,
        experiment_title = "Experiment Title",
        store_results = FALSE,
        quiet = FALSE,
        parallel = FALSE,
        BPPARAM = BiocParallel::bpparam()){
    
    if(!quiet) message("Initiating DEXSeq")
    
    stopifnot(
        "Please ensure the rownames have the format 'gene:transcript'" =
            grepl(":", rownames(dds)[1]))
    
    # Extract group and feature from rownames of DESeq2 object
    group_feat <- do.call("rbind", strsplit(rownames(dds), ":", fixed = TRUE))
    
    # Extract the found surrogate variables and covariates
    svs_covariates <- reduce_formula(
        DESeq2::design(dds), formularise = FALSE)
    
    # Add surrogate variables and covariates to DEXSeq design formula
    if(identical(svs_covariates, "1")){
        svs_covariates <- NULL
        design_formula <- formularise_vector(
            c("sample", "exon", "condition:exon"))
        reduced_formula <- formularise_vector(c("sample", "exon"))
    } else{
        design_formula <- formularise_vector(c(
            "sample", "exon", "condition:exon",
            paste0(svs_covariates, ":exon")))
        reduced_formula <- formularise_vector(c(
            "sample", "exon",
            paste0(svs_covariates, ":exon")))
        
    }
    
    # Convert to DEXSeqDataSet
    if(!quiet) message("Creating DEXSeqDataSet")
    sample_data <- create_sample_data(
        dds, group_col, baseline, case, svs_covariates)
    
    dds <- dds[, colnames(dds) %in% rownames(sample_data), drop = FALSE]
    
    dxd <- DEXSeq::DEXSeqDataSet(
        countData = DESeq2::counts(dds),
        sampleData = sample_data,
        design = design_formula,
        groupID = group_feat[, 1],
        featureID = group_feat[, 2]
    )
    
    # Store DEXSeqDataSet before DEXSeq analysis
    if(store_results) {
        store_result(
            dxd, paste0(experiment_title, "_dxd.RDS"),
            "DEXSeqDataSet", quiet = quiet)
    }
    
    ### Run DEXSeq
    if(!quiet) message("\nRunning DEXSeq -- This might take a while")
    if(!parallel) BiocParallel::register(BiocParallel::SerialParam())
    splicing_results <- DEXSeq::DEXSeq(
        dxd,
        reducedModel = reduced_formula,
        BPPARAM = BPPARAM,
        quiet = quiet)
    
    # Store result extraction
    if(store_results) {
        store_result(
            splicing_results,
            paste0(experiment_title, "_splicing_results.RDS"),
            "differential splicing results", quiet = quiet)
    }
    
    # Rename LFC column for consistency and human-readability
    # splicing_results <- data.frame(splicing_results)
    splicing_results <- splicing_results[, c(
        "groupID", "featureID", "exonBaseMean",
        "log2fold_C_B", "padj", "pvalue")]

    # Rename columns for consistency and human-readability
    colnames(splicing_results) <- c(
        "gene", "transcript", "baseMean", "lfc", "padj", "pvalue")
    
    return(splicing_results)
}


#' Run limma and voom analyses
#' 
#' This internal function runs a differential gene expression analysis using
#' limma (See \code{\link[limma]{lmFit}} for more detailed information), and 
#' differential gene splicing using voom (See \code{\link[limma]{voom}}).
#' Here, the surrogate variables found by \code{\link{run_sva}},
#' if any, will be added to the DESeqDataSet before running the analysis.
#'   
#' @inheritParams paired_diff
#' @inheritParams run_sva
#' @importFrom limma voom lmFit eBayes topTable
#' diffSplice topSplice is.fullrank
#' @import DESeq2
#' @importFrom stats model.matrix
#' @importFrom SummarizedExperiment colData
#' @keywords internal
#' @noRd
#' @return A list of two DFrames
#' @usage
#' run_limma(
#'     dds,
#'     group_col,
#'     baseline,
#'     case,
#'     experiment_title = "Experiment Title",
#'     store_results = FALSE,
#'     quiet = FALSE
#'     )
#' 
run_limma <- function(
        dds,
        group_col,
        baseline,
        case,
        experiment_title = "Experiment Title",
        store_results = FALSE,
        quiet = FALSE
){
    
    if(!quiet) message("Running limma+voom")
    
    # Design matrix
    modelMatrix <- stats::model.matrix(
        DESeq2::design(dds),
        data = SummarizedExperiment::colData(dds))
    
    stopifnot(
        "Experiment design is not full rank" 
        = limma::is.fullrank( modelMatrix ))
    
    # Voom
    voom <- limma::voom(
        counts = DESeq2::counts(dds),
        design = modelMatrix
    )
    fit <- limma::lmFit(object = voom)
    
    # Differential expression
    dgeModel <- limma::eBayes(fit)
    expression <- limma::topTable(
        dgeModel, sort.by = "p", coef = 2, number = Inf)
    
    # Differential splicing
    dguModel <- limma::diffSplice( 
        fit,
        exonid = sub(".*:", "", rownames(fit)),
        geneid = sub(":.*", "", rownames(fit)),
        verbose = !quiet
    )
    splicing <- limma::topSplice(
        dguModel, coef = 2, number = Inf, sort.by = "p", test = "t")
    
    
    # Convert result to data.frame and separate gene and transcript
    # expression <- data.frame(expression, row.names = rownames(expression))
    expression$gene <- sub(":.*", "", rownames(expression))
    expression$transcript <- sub(".*:", "", rownames(expression))
    names(expression) <- c(
        "lfc", "baseMean", "t", "pvalue", "padj", "B", "gene", "transcript")
    expression <- expression[, c(
        "gene", "transcript", "lfc", "padj", "pvalue", "baseMean")]
    
    # splicing <- data.frame(splicing, row.names = NULL)
    names(splicing) <- c("transcript", "gene", "lfc", "Test", "pvalue", "padj")
    splicing <- splicing[, c("gene", "transcript", "lfc", "padj", "pvalue")]
    
    
    baseMean <- expression$baseMean
    names(baseMean) <- expression$transcript
    splicing$baseMean <- baseMean[as.character(splicing$transcript)]
    
    # Store results before aggregation
    if(store_results) {
        store_result(
            list("expression" = expression, "splicing" = splicing),
            paste0(experiment_title, "_limma_results.RDS"),
            "differential expression/splicing results", quiet = quiet)
    }
    
    return(list(
        "expression" = S4Vectors::DataFrame(expression),
        "splicing" = S4Vectors::DataFrame(splicing)))
}



#' p-value aggregation to gene level
#' 
#' This internal function aggregates p-values to gene level using
#' Lancaster aggregation.
#' 
#' @param df A data.frame-type object with pvalues to aggregate
#' @param gene (Default: "gene") The column with gene names
#' @param p (Default: "pvalue") The column with p-values
#' @param lfc (Default: "lfc") The column with log2-foldchanges
#' @param type (Default: "expression") 
#' Either "expression" for aggregation of
#' DESeq2 results or "splicing" for aggregation of DEXSeq results
#' @param weights (Default: "baseMean") The column to use
#' for weighting the aggregation
#' @importFrom aggregation lancaster
#' @importFrom stats weighted.mean
#' @importFrom S4Vectors split DataFrame
#' @keywords internal
#' @noRd
#' @usage 
#' aggregate_pvalue(
#'     df,
#'     type = c("expression", "splicing")
#'     )
#' @return A DFrame
aggregate_pvalue <- function(
        df,
        type = c("expression", "splicing")) {
    
    # Check input data
    required_cols <- c("gene", "pvalue", "baseMean", "lfc")
    stopifnot(
        "Input data is not a data.frame" = 
            is(df, "data.frame") | is(df, "DFrame"),
        "Columns 'gene', 'pvalue', 'baseMean', and 'lfc' are not in input data"
        = all(required_cols %in% colnames(df)))
    
    # Check input parameters
    type <- match.arg(type)
    
    # Rename columns and convert values to the correct type
    df <- df[!is.na(df$pvalue) & !is.na(df$baseMean),]
    # Prevent warning from Lancaster
    if(any(df$pvalue < 10e-320)) {
        if (is(df, "data.frame")) df$pvalue[df$pvalue < 10e-320] <- 10e-320
        else df@listData$pvalue[df@listData$pvalue < 10e-320] <- 10e-320
        }

    # Split data by gene and perform aggregation
    res <- S4Vectors::split(df, df$gene)
    if(type == "splicing") {
        res <- lapply(res, function(x) {
            lfc <- x$lfc[which.min(x$pvalue)]
            pvalue <- aggregation::lancaster(x$pvalue, weights = x$baseMean)
            baseMean <- sum(x$baseMean)
            data.frame(lfc = lfc, pvalue = pvalue, baseMean = baseMean)
        })
    } else if(type == "expression") {
        res <- lapply(res, function(x) {
            lfc <- weighted.mean(x$lfc, w = x$baseMean)
            pvalue <- aggregation::lancaster(x$pvalue, weights = x$baseMean)
            baseMean <- sum(x$baseMean)
            data.frame(lfc = lfc, pvalue = pvalue, baseMean = baseMean)
        })
    }
    
    # Combine results and adjust p values
    res <- do.call("rbind", res)
    res <- S4Vectors::DataFrame(gene = rownames(res), res)
    res$padj <- stats::p.adjust(res$pvalue, "fdr")
    
    return(res)
}
