data("example_diff_result")
data("example_gene_sets")

example_diff_result_with_dif <- example_diff_result
example_diff_result_with_dif$lfc_expression <- 1
example_diff_result_with_dif$max_abs_dif_splicing <- 1



test_that("paired_ora works", {
    ora_test <- paired_ora(example_diff_result_with_dif, example_gene_sets)
    expect_equal(nrow(ora_test), 8)
    ora_test <- paired_ora(
        example_diff_result, example_gene_sets, expression_only = TRUE)
    expect_gt(nrow(ora_test), 0)
    ora_test <- paired_ora(
        example_diff_result, example_gene_sets, effect_size_filter = FALSE)
    expect_equal(nrow(ora_test), 8)
})

test_that("paired_ora prints results", {
    ora_test <- paired_ora(
        example_diff_result_with_dif,
        example_gene_sets,
        experiment_title = "test")
    expect_true(file.exists("results/test_ora.RDS"))
    file.remove("results/test_ora.RDS")
    expect_false(file.exists("results/test_ora.RDS"))
    
    ora_test <- paired_ora(example_diff_result, example_gene_sets, experiment_title = "test", expression_only = TRUE)
    expect_true(file.exists("results/test_ora.RDS"))
    file.remove("results/test_ora.RDS")
    expect_false(file.exists("results/test_ora.RDS"))
})

test_that("paired_ora works with limma", {
    data("example_se")
    diff_results <- suppressWarnings(
        paired_diff(
            object = example_se,
            group_col = "group_nr",
            sample_col = "id",
            baseline = "1",
            case = "2",
            use_limma = TRUE
        ))
    
    ora_test <- paired_ora(
        diff_results, example_gene_sets, effect_size_filter = FALSE)
    expect_equal(nrow(ora_test), 8)
})

test_that("prepare_msigdb returns a list of gene sets", {
    gene_sets <- prepare_msigdb()
    expect_true(is.list(gene_sets))
    expect_true(length(gene_sets) > 0)
    expect_true(all(sapply(gene_sets, is.character)))
})

test_that("it subsets genes to a cutoff based on type", {
    expression_genes <- subset_genes(example_diff_result, "expression", 0.05)
    expect_true(all(expression_genes$padj_expression < 0.05))
    
    splicing_genes <- subset_genes(example_diff_result, "splicing", 0.05)
    expect_true(all(splicing_genes$padj_splicing < 0.05))
})

test_that("effect-size filtering uses expression LFC and splicing dIF", {
    toy_diff <- data.frame(
        gene = paste0("g", 1:4),
        pvalue_expression = rep(0.01, 4),
        padj_expression = c(0.01, 0.01, 0.2, 0.2),
        lfc_expression = c(0.1, 1, 1, 1),
        pvalue_splicing = rep(0.01, 4),
        padj_splicing = c(0.2, 0.2, 0.01, 0.01),
        max_abs_dif_splicing = c(0.2, 0.2, 0.05, 0.2)
    )
    
    expression_genes <- subset_genes(
        toy_diff, "expression", 0.05, effect_size_filter = TRUE)
    expect_equal(expression_genes$gene, "g2")
    
    splicing_genes <- subset_genes(
        toy_diff, "splicing", 0.05, effect_size_filter = TRUE)
    expect_equal(splicing_genes$gene, "g4")
    
    paired_genes <- subset_genes(
        toy_diff, "paired", 0.05, effect_size_filter = TRUE)
    expect_equal(paired_genes$gene, c("g2", "g4"))
})

test_that("effect-size filtering does not change ORA universe", {
    toy_diff <- data.frame(
        gene = paste0("g", 1:4),
        pvalue_expression = rep(0.01, 4),
        padj_expression = c(0.01, 0.01, 0.2, 0.2),
        lfc_expression = c(0.1, 1, 1, 1),
        pvalue_splicing = rep(0.2, 4),
        padj_splicing = rep(0.2, 4),
        max_abs_dif_splicing = rep(0, 4)
    )
    gene_sets <- list(target = c("g1", "g2"))
    
    ora_test <- run_ora(
        toy_diff,
        gene_sets = gene_sets,
        type = "expression",
        cutoff = 0.05,
        min_size = 1,
        effect_size_filter = TRUE)
    
    expect_equal(ora_test$overlap, 1)
    expect_equal(ora_test$size, 2)
    expect_equal(ora_test$relative_risk, 2)
})

test_that("missing splicing dIF gives a useful error", {
    diff_without_dif <- as.data.frame(example_diff_result)
    diff_without_dif <- diff_without_dif[
        , colnames(diff_without_dif) != "max_abs_dif_splicing"]
    
    expect_error(
        paired_ora(diff_without_dif, example_gene_sets),
        "max_abs_dif_splicing")
})

# Create a mock ORA data frame
ora <- data.frame(
    pathway = c("Pathway A", "Pathway B", "Pathway C"),
    overlap = c(10, 10, 10),
    size = c(100, 200, 300),
    stringsAsFactors = FALSE
)

# Define test cases
test_that("compute_enrichment returns correct enrichment scores", {
    # Test with n_genes = 1000 and n_universe = 10000
    result <- compute_enrichment(ora, n_genes = 1000, n_universe = 10000)
    expect_equal(result$relative_risk, c(1, 0.5, 1/3))
    expect_equal(result$enrichment_score, c(log2(1+0.06), log2(0.5+0.06), log2(1/3+0.06)))
})
