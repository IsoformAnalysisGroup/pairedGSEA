#' Wrapper for missing packages
#'
#' @noRd
check_missing_package <- function(package, repo = "CRAN", git_repo = ""){
    
    if (repo == "CRAN"){
        install_function <- "install.packages('"
    } else if (repo == "github") {
        install_function <- paste0("devtools::install_github('", git_repo, "/")
    } else if (repo == "Bioc"){
        install_function <- "BiocManager::install('"
    } else{
        install_function <- "Unknown repository.. "
    }
    
    if(!requireNamespace(package, quietly = TRUE)){
        stop(
            "Package ", package," is not installed.\n",
            "Please run: ", install_function, package, "')")
    }
    requireNamespace(package, quietly = TRUE)
}

#' Check if directory exists, if not, make it
#' @noRd
check_make_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {dir.create(dir_path)}
}

#' Check colname
#' @noRd
check_colname <- function(df_colnames, col_name, location = "metadata"){
    if(is(df_colnames, "data.frame") | is(df_colnames, "DFrame")) 
        df_colnames <- colnames(df_colnames)
    if(!is.null(col_name)){
        if(!(col_name %in% df_colnames)){
            stop("Column \"", col_name, "\" was not found in the ", location)
    }}
}

#' Check comparison has the right format
#' @noRd
check_comparison <- function(comparison){
    
    if(is(comparison, "list")) comparison <- as.character(comparison)
    if(is(comparison, "character") & length(comparison) == 1){
        comparison <- gsub("vs", "v", comparison)
        stopifnot(
            "Comparison must have the format 'baseline_v_case' (e.g., '2v1')." =
                grepl("v", comparison))
        stopifnot(
            "Comparison must not contain multiple 'v's, the format should be
            'baseline_v_case' (e.g., '2v1')." =
                length(gregexpr("v", comparison)[[1]]) == 1)
        comparison <- gsub(pattern = " ", replacement = "", comparison)
        comparison <- gsub(pattern = "_v_", replacement = "v", comparison)
        comparison <- strsplit(comparison, "v")[[1]]
    }
    stopifnot(
        "Comparison must be a character string (e.g., '2v1')
        or a list (e.g., c('2', '1')). The format being 'baseline_v_case'." =
            is(comparison, "character") & (length(comparison) == 2))
    return(comparison)
}


#' Store result object
#' @importFrom utils write.csv write.table
#' @note Suggested: importFrom writexl write_xlsx
#' @noRd
store_result <- function(
        object, file, analysis = "results", quiet = FALSE,
        location = "results/"){
    check_make_dir(location)
    if(!startsWith(file, location)) file <- paste0(location, file)
    
    if(endsWith(toupper(file), ".RDS")) saveRDS(object, file = file)
    else if(endsWith(tolower(file), ".rdata")) save(object, file = file)
    else if(endsWith(tolower(file), ".csv")) utils::write.csv(
        object, file = file)
    else if(endsWith(tolower(file), ".tsv")) utils::write.table(
        object, file = file)
    else if(endsWith(tolower(file), ".xlsx")) {check_missing_package("writexl")
        writexl::write_xlsx(object, path = file)}
    else{stop("Could not store result object at ", file)}
    if(!quiet) message("Stored ", analysis, " in ", file)
}



#' Pre-filter
#' @inheritParams paired_diff
#' @noRd
pre_filter <- function(dds, threshold = 10, quiet = FALSE){
    if(threshold < 1) return(dds)
    # Remove low counts
    remove_low <- rowSums(DESeq2::counts(dds)) < threshold
    dds <- dds[!remove_low,]
    
    # Remove rows with counts in only one sample
    remove_singles <- rowSums(DESeq2::counts(dds) > 0) < 2
    
    dds <- dds[!remove_singles,]
    
    if(sum(remove_low) > 0 | sum(remove_singles) > 0 & !quiet) message(
        "\nRemoving ", sum(remove_low), 
        " rows with a summed count lower than ", threshold,
        "\nRemoving ", sum(remove_singles),
        " rows with counts in less than 2 samples.\n")
    
    return(dds)
}

#' Convert character vector into a stats formula
#' @param vector Character vector of variables to add to model
#' @importFrom stats formula
#' @noRd
formularise_vector <- function(vector){
    stopifnot(
        "Names must not contain spaces" = all(grepl("\\s", vector) == FALSE))
    if(is.null(vector)){
        formula <- NULL
    } else if(length(vector) == 1){
        formula <- stats::formula(paste0("~", vector))
    } else{
        formula <- stats::formula(paste0(
            "~", vector[1], 
            paste("+", vector[2:length(vector)], collapse = " ")))
    }
    return(formula)
}

#' Reduce a design formula by removing one variable
#' 
#' @param formula Design formula to reduce.
#' @param formularise (Default: TRUE) Logical determining if the design
#' should be formularised
#' @param term Optional term to remove. If omitted, the first term is removed
#' for backwards compatibility.
#' @noRd
reduce_formula <- function(formula, formularise = TRUE, term = NULL){
    formula_vector <- attr(stats::terms(formula), "term.labels")

    if(length(formula_vector) == 0) {
        reduced_vector <- "1"
    } else{
        if(is.null(term)) term <- formula_vector[1]
        if(length(term) != 1 || !(term %in% formula_vector)) {
            stop("Term '", term, "' was not found in the design formula.",
                 call. = FALSE)
        }
        reduced_vector <- formula_vector[formula_vector != term]
        if(length(reduced_vector) == 0) reduced_vector <- "1"
    }
    
    if(formularise) return(formularise_vector(reduced_vector))
    
    return(reduced_vector)
}

#' Create sample data for differential expression analysis
#'
#' @param dds A DESeq2 object.
#' @param group_col A character string specifying the column in the
#'  colData to use for grouping samples.
#' @param baseline A character string specifying the value in \code{group_col}
#'  to use as the baseline group.
#' @param case A character string specifying the value in \code{group_col}
#'  to use as the case group.
#' @param svs_covariates A character vector specifying the columns to
#'  include in the output data.frame in addition to the \code{group_col} column.
#' @noRd
#' @keywords internal
#' @return A data.frame with the formatted sample data.
#'
#' @examples
#' create_sample_data(dds, group_col = "condition", baseline = "1", case = "2",
#'  svs_covariates = c("sv1", "sv2"))
#'
create_sample_data <- function(
        dds, group_col, baseline, case, svs_covariates) {
    # Extract colData from DESeq2 object and convert to data.frame
    sample_data <- as.data.frame(SummarizedExperiment::colData(dds))
    
    # Set values in group_col column to B or C based on whether they match 
    # the baseline or case variables
    sample_data[[group_col]] <- ifelse(
        sample_data[[group_col]] == baseline, "B", 
        ifelse(sample_data[[group_col]] == case, "C", NA))
    
    # Select columns specified in svs_covariates
    sample_data <- sample_data[, c(group_col, svs_covariates), drop = FALSE]
    
    # Remove any rows with missing values and convert group_col
    #  to a factor with levels B and C
    sample_data <- na.omit(sample_data)
    sample_data[[group_col]] <- factor(
        sample_data[[group_col]], levels = c("B", "C"))
    
    # Rename group_col to "condition"
    colnames(sample_data)[which(
        colnames(sample_data) == group_col)] <- "condition"
    
    return(sample_data)
}
