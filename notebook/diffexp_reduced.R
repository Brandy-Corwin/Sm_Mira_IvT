library(tidyverse)
library(broom)
library(DESeq2)





# get rid of age and duplicate for hisat
# change more as needed




star_counts_df <- read_tsv("counting_reduced/counts/star_counts.tsv", comment = "#") |>
             mutate(across(where(is.numeric), as.integer))

star_counts_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/star.bam:"), everything()) |>
    rowwise() |>
    mutate(total_counts = sum(c_across(where(is.numeric)), na.rm = T)) |>
    filter(total_counts >= 10)


star_sample_summary <- star_counts_df |>
    select(Geneid, contains('star.bam')) |>
    rename_with(~str_remove(., "counting_reduced/dedup/star.bam:"), everything()) |>
    pivot_longer(-Geneid, names_to = 'sample', values_to = 'count') |>
    filter(count > 0) |>
    group_by(Geneid) |>
    tally() |>
    filter(n <= 3)


star_genes_to_remove = star_sample_summary$Geneid

star_counts_filt <- star_counts_summary |>
    filter(!Geneid %in% genes_to_remove) |>
    arrange(Geneid) |>
    select(-total_counts)


star_counts_m <- star_counts_filt |>
    select(-Geneid) |>
    as.matrix()
rownames(star_counts_m) <- star_counts_filt$Geneid


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


star_pca_fit <- t(log10(star_counts_m + 1)) |> 
  prcomp(scale = TRUE)


star_pca_fit |>
  augment(t(star_counts_m)) |>
  dplyr::rename(sample = .rownames) |>
  separate(sample, into = c('tissue', 'age'), sep = '_') |>
  mutate(age = str_remove(age, '[0-9]')) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = tissue, shape = age)) + 
  geom_point(size = 4)

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
star_volcano_plot
