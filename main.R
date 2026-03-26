library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts_df <- read.delim(counts_csv, stringsAsFactors = FALSE, check.names = FALSE)
  meta_df <- read.csv(metafile_csv, stringsAsFactors = FALSE)
  
  meta_df <- meta_df %>%
    dplyr::filter(timepoint %in% selected_times)
  
  meta_df$timepoint <- factor(meta_df$timepoint, levels = selected_times)
  meta_df$timepoint <- relevel(meta_df$timepoint, ref = "vP0")
  
  counts_mat <- counts_df[, meta_df$samplename, drop = FALSE]
  counts_mat <- as.matrix(counts_mat)
  rownames(counts_mat) <- counts_df[[1]]
  
  coldata <- meta_df %>%
    dplyr::select(samplename, timepoint)
  
  rownames(coldata) <- coldata$samplename
  
  se <- SummarizedExperiment(
    assays = list(counts = counts_mat),
    colData = coldata
  )
  
  return(se)
}

return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res_df <- as.data.frame(results(dds))
  
  return(list(res_df, dds))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  labeled_results <- deseq2_res %>%
    as_tibble(rownames = "genes") %>%
    mutate(
      volc_plot_status = case_when(
        !is.na(padj) & padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        !is.na(padj) & padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  
  return(labeled_results)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(bins = 30) +
    theme_minimal() +
    labs(
      x = "pvalue",
      y = "count",
      title = "Histogram of DESeq2 p-values"
    )
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  plot_df <- labeled_results %>%
    filter(!is.na(padj), padj < padj_threshold)
  
  ggplot(plot_df, aes(x = log2FoldChange)) +
    geom_histogram(bins = 30) +
    theme_minimal() +
    labs(
      x = "log2FoldChange",
      y = "count",
      title = "Histogram of significant log2 fold changes"
    )
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- labeled_results %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = num_genes) %>%
    pull(genes)
  
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  plot_df <- as.data.frame(norm_counts[top_genes, , drop = FALSE]) %>%
    tibble::rownames_to_column(var = "genes") %>%
    pivot_longer(
      cols = -genes,
      names_to = "samplename",
      values_to = "normalized_count"
    )
  
  meta_df <- as.data.frame(colData(dds_obj)) %>%
    tibble::rownames_to_column(var = "samplename")
  
  plot_df <- left_join(plot_df, meta_df, by = "samplename")
  
  ggplot(plot_df, aes(x = timepoint, y = normalized_count)) +
    geom_point(position = position_jitter(width = 0.1, height = 0)) +
    facet_wrap(~ genes, scales = "free_y") +
    theme_minimal() +
    labs(
      x = "timepoint",
      y = "normalized count",
      title = "Normalized counts for top genes"
    )
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  plot_df <- labeled_results %>%
    mutate(neg_log10_padj = -log10(padj))
  
  ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = volc_plot_status)) +
    geom_point() +
    theme_minimal() +
    labs(
      x = "log2FoldChange",
      y = "-log10(padj)",
      color = "status",
      title = "Volcano plot"
    )
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read.delim(id2gene_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(id2gene) <- c("genes", "symbol")
  
  ranked_df <- labeled_results %>%
    dplyr::select(genes, log2FoldChange) %>%
    left_join(id2gene, by = "genes") %>%
    filter(!is.na(symbol), !is.na(log2FoldChange)) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    arrange(desc(log2FoldChange))
  
  rnk_list <- ranked_df$log2FoldChange
  names(rnk_list) <- ranked_df$symbol
  
  return(rnk_list)
}
#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  
  fgsea_results <- fgsea(
    pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  ) %>%
    as_tibble()
  
  return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top_pos <- fgsea_results %>%
    arrange(desc(NES)) %>%
    slice_head(n = num_paths)
  
  top_neg <- fgsea_results %>%
    arrange(NES) %>%
    slice_head(n = num_paths)
  
  plot_df <- bind_rows(top_pos, top_neg) %>%
    mutate(pathway = forcats::fct_reorder(pathway, NES))
  
  ggplot(plot_df, aes(x = pathway, y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(
      x = "pathway",
      y = "NES",
      title = "Top enriched pathways"
    )
}

