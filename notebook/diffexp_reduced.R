library(tidyverse)
library(broom)
library(DESeq2)



# get rid of age


star_counts_df <- read_tsv("counting_reduced/counts/star_counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
write.csv(star_counts_df, "count_matrices_reduced/star_counts_df.csv")

star_counts_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
write.csv(star_counts_summary, "count_matrices_reduced/star_counts_summary.csv")

star_sample_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)
write.csv(star_sample_summary, "count_matrices_reduced/star_sample_summary.csv")

star_genes_to_remove = star_sample_summary$Geneid

star_counts_filt <- star_counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

star_counts_m <- star_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(star_counts_m) <- star_counts_filt$Geneid
write.csv(star_counts_m, "count_matrices_reduced/star_counts_m.csv")

star_dists <- dist(t(star_counts_m))

star_dists_df <- as.matrix(star_dists) |>
    as_tibble(rownames = 'sample')

star_dist_plot <- star_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
ggsave("plots_reduced/star_dist_plot.png")

star_pca_fit <- t(log10(star_counts_m + 1)) |> 
  prcomp(scale = TRUE)


# age starts here: fix below
star_pca_fit |>
  augment(t(star_counts_m)) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'age'), sep = '_') |>
  mutate(age = str_remove(age, '[0-9]')) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = tissue, shape = age)) + 
  geom_point(size = 4)
ggsave("plots_reduced/star_pca_fit.png")

star_metadata <- data.frame(sample_id = colnames(star_counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           age = str_sub(sample_id, 5, 6),
           rep = str_sub(sample_id, 7))
rownames(star_metadata) <- star_metadata$sample_id
star_metadata <- select(star_metadata, -sample_id)


all(rownames(star_metadata) == colnames(star_counts_m))

star_metadata2 <- star_metadata |>
    mutate(tissue_age = str_c(tissue, "_", age)) |>
    select(tissue_age, rep)


star_dds2 <- DESeqDataSetFromMatrix(countData = star_counts_m,
                               colData = star_metadata2,
                               design = ~ tissue_age)
star_dds2 <- star_DESeq(dds2)


star_LIV_vs_INT <- results(star_dds2, contrast = c('tissue_age', 'LIV_ma', 'INT_ma'))


as_tibble(star_LIV_vs_INT, rownames = 'gene_id') |>
    filter(gene_id == 'Smp_333930')

star_volcano_data <- as_tibble(star_LIV_vs_INT, rownames = 'gene_id')

star_volcano_plot <- star_volcano_data |> 
    ggplot() +
    geom_point(aes(x = log2FoldChange, y = -log10(padj)), 
               color = '#CD5C5C', size = 5) +
    geom_hline(yintercept = -log10(0.05))
ggsave("plots_reduced/star_volcano_plot.png")






hisat_counts_df <- read_tsv("counting_reduced/counts/hisat_counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))
write.csv(hisat_counts_df, "count_matrices_reduced/hisat_counts_df.csv")

hisat_counts_summary <- hisat_counts_df |>
    select(Geneid, contains('hisat.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/hisat.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)
write.csv(hisat_counts_summary, "count_matrices_reduced/hisat_counts_summary.csv")

hisat_sample_summary <- hisat_counts_df |>
    select(Geneid, contains('hisat.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/hisat.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)
write.csv(hisat_sample_summary, "count_matrices_reduced/hisat_sample_summary.csv")

hisat_genes_to_remove = hisat_sample_summary$Geneid

hisat_counts_filt <- hisat_counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)

hisat_counts_m <- hisat_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(hisat_counts_m) <- hisat_counts_filt$Geneid
write.csv(hisat_counts_m, "count_matrices_reduced/hisat_counts_m.csv")

hisat_dists <- dist(t(hisat_counts_m))

hisat_dists_df <- as.matrix(hisat_dists) |>
    as_tibble(rownames = 'sample')

hisat_dist_plot <- hisat_dists_df |>
    pivot_longer(-sample, names_to = 'comp', values_to = 'dist') |>
    ggplot(aes(x = sample, y = comp, fill = dist)) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_equal() +
    NULL
ggsave("plots_reduced/hisat_dist_plot.png")

hisat_pca_fit <- t(log10(hisat_counts_m + 1)) |> 
  prcomp(scale = TRUE)


# age starts here: fix below
hisat_pca_fit |>
  augment(t(hisat_counts_m)) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'age'), sep = '_') |>
  mutate(age = str_remove(age, '[0-9]')) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = tissue, shape = age)) + 
  geom_point(size = 4)
ggsave("plots_reduced/hisat_pca_fit.png")

hisat_metadata <- data.frame(sample_id = colnames(hisat_counts_m)) |>
    mutate(tissue = str_sub(sample_id, 1, 3),
           age = str_sub(sample_id, 5, 6),
           rep = str_sub(sample_id, 7))
rownames(hisat_metadata) <- hisat_metadata$sample_id
hisat_metadata <- select(hisat_metadata, -sample_id)


all(rownames(hisat_metadata) == colnames(hisat_counts_m))

hisat_metadata2 <- hisat_metadata |>
    mutate(tissue_age = str_c(tissue, "_", age)) |>
    select(tissue_age, rep)


hisat_dds2 <- DESeqDataSetFromMatrix(countData = hisat_counts_m,
                               colData = hisat_metadata2,
                               design = ~ tissue_age)
hisat_dds2 <- hisat_DESeq(dds2)


hisat_LIV_vs_INT <- results(hisat_dds2, contrast = c('tissue_age', 'LIV_ma', 'INT_ma'))


as_tibble(hisat_LIV_vs_INT, rownames = 'gene_id') |>
    filter(gene_id == 'Smp_333930')

hisat_volcano_data <- as_tibble(hisat_LIV_vs_INT, rownames = 'gene_id')

hisat_volcano_plot <- hisat_volcano_data |> 
    ggplot() +
    geom_point(aes(x = log2FoldChange, y = -log10(padj)), 
               color = '#CD5C5C', size = 5) +
    geom_hline(yintercept = -log10(0.05))
ggsave("plots_reduced/hisat_volcano_plot.png")