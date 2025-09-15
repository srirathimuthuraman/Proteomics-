
#### Male Proteomics Analysis ####
library(openxlsx)
library(proDA)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(vegan)

# Load data
filtered_data_M <- as.matrix(read.xlsx("Filtered_data_M_gene_names_Final.xlsx", rowNames = T))

# Filter proteins present in at least 6 samples per group
presence_threshold <- 6
MKO_samples <- grep("MKO_", colnames(filtered_data_M), value = TRUE)
MWT_samples <- grep("MWT_", colnames(filtered_data_M), value = TRUE)

MKO_present <- rowSums(!is.na(filtered_data_M[, MKO_samples]) & filtered_data_M[, MKO_samples] != 0) >= presence_threshold
MWT_present <- rowSums(!is.na(filtered_data_M[, MWT_samples]) & filtered_data_M[, MWT_samples] != 0) >= presence_threshold

filtered_data_M1 <- filtered_data_M[MKO_present | MWT_present, ]

# Automatically assign group labels
groups_M <- factor(gsub("_.*", "", colnames(filtered_data_M1)))

#SETSEED
set.seed(13) 

# Imputation function
impute_group <- function(group_row) {
  imputed <- unlist(group_row)
  non_na_vals <- imputed[!is.na(imputed) & imputed != 0]
  if (length(non_na_vals) > 0) {
    min_val <- min(non_na_vals)
    imputed[is.na(imputed)] <- runif(sum(is.na(imputed)), min = min_val * 0.9, max = min_val * 1.1)
  }
  return(imputed)
}

# Apply imputation
data_processed <- filtered_data_M1
for (i in 1:nrow(data_processed)) {
  data_processed[i, MKO_samples] <- impute_group(data_processed[i, MKO_samples])
  data_processed[i, MWT_samples] <- impute_group(data_processed[i, MWT_samples])
}

# Ensure no NAs
data_processed[is.na(data_processed)] <- 0
data_processed_M <- data_processed

# Differential analysis using proDA
fit_M <- proDA(data_processed_M, design = groups_M)
prot_diff_M_out <- proDA::test_diff(fit_M, MKO - MWT, pval_adjust_method = "BH")

# Add -log10 p-value and gene regulation category
prot_diff_M_out$log_pval <- -log10(prot_diff_M_out$pval)
prot_diff_M_out <- prot_diff_M_out[prot_diff_M_out$name != "Negr1", ]
prot_diff_M_out$regulation <- with(prot_diff_M_out,
                                   ifelse(pval < 0.05 & diff > 0, "Upregulated",
                                          ifelse(pval < 0.05 & diff < 0, "Downregulated", "Not Significant")))

# Count regulated proteins
upregulated_count <- sum(prot_diff_M_out$regulation == "Upregulated")
downregulated_count <- sum(prot_diff_M_out$regulation == "Downregulated")

# Top 10 regulated proteins (only if observed in all 16 samples)
filtered_proteins <- prot_diff_M_out %>% filter(n_obs == 16)
top_upregulated <- filtered_proteins %>% filter(regulation == "Upregulated") %>% arrange(desc(diff)) %>% slice_head(n = 10)
top_downregulated <- filtered_proteins %>% filter(regulation == "Downregulated") %>% arrange(diff) %>% slice_head(n = 10)
top_proteins <- bind_rows(top_upregulated, top_downregulated)

# To explort table 
openxlsx::write.xlsx(prot_diff_M_out, "prot_diff_M_out.xlsx")

# Load data
file_path <- "Prot_diff_M_out.xlsx"
prot_diff_M_out <- read.xlsx(file_path, rowNames = TRUE)

# Remove "Negr1" and handle NAs
prot_diff_M_out <- prot_diff_M_out[!rownames(prot_diff_M_out) %in% "Negr1", ]
prot_diff_M_out <- prot_diff_M_out %>% filter(!is.na(pval) & !is.na(diff))

# Convert types
prot_diff_M_out$pval <- as.numeric(prot_diff_M_out$pval)
prot_diff_M_out$diff <- as.numeric(prot_diff_M_out$diff)

# Add log p-value and regulation status
prot_diff_M_out$log_pval <- -log10(prot_diff_M_out$pval)
prot_diff_M_out$regulation <- with(prot_diff_M_out,
                                   ifelse(pval < 0.05 & diff > 0, "Upregulated",
                                          ifelse(pval < 0.05 & diff < 0, "Downregulated", "Not Significant")))

# Count regulation
upregulated_count <- sum(prot_diff_M_out$regulation == "Upregulated")
downregulated_count <- sum(prot_diff_M_out$regulation == "Downregulated")

# Top 12 proteins by p-value
top_12_upregulated <- prot_diff_M_out %>%
  filter(regulation == "Upregulated") %>%
  arrange(pval) %>%
  slice_head(n = 12)

top_12_downregulated <- prot_diff_M_out %>%
  filter(regulation == "Downregulated") %>%
  arrange(pval) %>%
  slice_head(n = 12)


# Plot
ggplot(prot_diff_M_out, aes(x = diff, y = log_pval)) +
  geom_point(aes(color = regulation), size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  labs(
    title = paste0("Proteins (", upregulated_count, " Upregulated, ", downregulated_count, " Downregulated)"),
    x = "Log2 Fold Change (Difference in Abundance)",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(-2.0, 2.5) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),  # Larger plot title
    axis.title.x = element_text(size = 14),             # Increase x-axis label size
    axis.title.y = element_text(size = 14),             # Increase y-axis label size
    axis.text = element_text(size = 12)                 # Increase axis tick text size
  ) +
  geom_text_repel(
    data = top_12_upregulated,
    aes(label = rownames(top_12_upregulated)),
    size = 6,
    box.padding = 0.4,
    max.overlaps = 100
  ) +
  geom_text_repel(
    data = top_12_downregulated,
    aes(label = rownames(top_12_downregulated)),
    size = 6,
    box.padding = 0.4,
    max.overlaps = 100
  )

# Step 1: Filter out "Negr1"
data_processed_M_out <- data_processed_M[rownames(data_processed_M) != "Negr1", ]

# Step 2: Log2 transform the data (after adding 1 to avoid -Inf)
data_pca <- t(log2(data_processed_M_out + 1))

# Step 3: Run PCA
pca_res <- prcomp(data_pca, scale. = FALSE)

# Step 4: Calculate percent variance explained
var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# Step 5: Create PCA dataframe
pca_df <- as.data.frame(pca_res$x)

# Extract group info from row names (e.g., "FKO_1", "FWT_2", etc.)
pca_df$Group <- factor(gsub("_.*", "", rownames(pca_df)))

# Step 6: Set custom colors
group_colors <- c("MKO" = "darkgreen", "MWT" = "purple")

# Step 7: Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 6, alpha = 0.8) +
  labs(
    title = "PCA of Male Protein Abundance",
    x = pc1_label,
    y = pc2_label
  ) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 12)
  )




#### Female Proteomics Analysis ####


# Load data
filtered_data_F <- as.matrix(read.xlsx("Filtered_data_F_gene_names_Final.xlsx", rowNames = T))
# Filter proteins present in at least 8 samples per group
presence_threshold <- 6
FKO_samples <- grep("FKO_", colnames(filtered_data_F), value = TRUE)
FWT_samples <- grep("FWT_", colnames(filtered_data_F), value = TRUE)

FKO_present <- rowSums(!is.na(filtered_data_F[, FKO_samples]) & filtered_data_F[, FKO_samples] != 0) >= presence_threshold
FWT_present <- rowSums(!is.na(filtered_data_F[, FWT_samples]) & filtered_data_F[, FWT_samples] != 0) >= presence_threshold

filtered_data_F1 <- filtered_data_F[FKO_present | FWT_present, ]

# Automatically assign group labels
groups_F <- factor(gsub("_.*", "", colnames(filtered_data_F1)))

#SETSEED
set.seed(13) 
# Imputation function
impute_group <- function(group_row) {
  imputed <- unlist(group_row)
  non_na_vals <- imputed[!is.na(imputed) & imputed != 0]
  if (length(non_na_vals) > 0) {
    min_val <- min(non_na_vals)
    imputed[is.na(imputed)] <- runif(sum(is.na(imputed)), min = min_val * 0.9, max = min_val * 1.1)
  }
  return(imputed)
}

# Apply imputation
data_processed <- filtered_data_F1
for (i in 1:nrow(data_processed)) {
  data_processed[i, FKO_samples] <- impute_group(data_processed[i, FKO_samples])
  data_processed[i, FWT_samples] <- impute_group(data_processed[i, FWT_samples])
}

# Ensure no NAs
data_processed[is.na(data_processed)] <- 0
data_processed_F <- data_processed

# Differential expression with proDA
fit_F <- proDA(data_processed_F, design = groups_F)
prot_diff_F_out <- proDA::test_diff(fit_F, FKO - FWT, pval_adjust_method = "BH")


# Add -log10(pval) and classification
prot_diff_F_out$log_pval <- -log10(prot_diff_F_out$pval)
prot_diff_F_out <- prot_diff_F_out[prot_diff_F_out$name != "Negr1", ]
prot_diff_F_out$regulation <- with(prot_diff_F_out,
                                   ifelse(pval < 0.05 & diff > 0, "Upregulated",
                                          ifelse(pval < 0.05 & diff < 0, "Downregulated", "Not Significant")))

# Count regulation types
upregulated_count <- sum(prot_diff_F_out$regulation == "Upregulated")
downregulated_count <- sum(prot_diff_F_out$regulation == "Downregulated")

# Top 10 up/downregulated proteins (in all samples)
filtered_proteins <- prot_diff_F_out %>% filter(n_obs == 16)
top_upregulated <- filtered_proteins %>% filter(regulation == "Upregulated") %>% arrange(desc(diff)) %>% slice_head(n = 10)
top_downregulated <- filtered_proteins %>% filter(regulation == "Downregulated") %>% arrange(diff) %>% slice_head(n = 10)
top_proteins <- bind_rows(top_upregulated, top_downregulated)

# To explort table 
openxlsx::write.xlsx(prot_diff_F_out, "prot_diff_F_out.xlsx")

# Load necessary libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Load data
file_path <- "Prot_diff_F_out.xlsx"
prot_diff_F_out <- read.xlsx(file_path, rowNames = TRUE)

# Remove "Negr1" and handle NAs
prot_diff_F_out <- prot_diff_F_out[!rownames(prot_diff_F_out) %in% "Negr1", ]
prot_diff_F_out <- prot_diff_F_out %>% filter(!is.na(pval) & !is.na(diff))

# Convert types
prot_diff_F_out$pval <- as.numeric(prot_diff_F_out$pval)
prot_diff_F_out$diff <- as.numeric(prot_diff_F_out$diff)

# Add log p-value and regulation status
prot_diff_F_out$log_pval <- -log10(prot_diff_F_out$pval)
prot_diff_F_out$regulation <- with(prot_diff_F_out,
                                   ifelse(pval < 0.05 & diff > 0, "Upregulated",
                                          ifelse(pval < 0.05 & diff < 0, "Downregulated", "Not Significant")))

# Count regulation
upregulated_count <- sum(prot_diff_F_out$regulation == "Upregulated")
downregulated_count <- sum(prot_diff_F_out$regulation == "Downregulated")

# Top 12 proteins by p-value
top_12_upregulated <- prot_diff_F_out %>%
  filter(regulation == "Upregulated") %>%
  arrange(pval) %>%
  slice_head(n = 12)

top_12_downregulated <- prot_diff_F_out %>%
  filter(regulation == "Downregulated") %>%
  arrange(pval) %>%
  slice_head(n = 12)


# Plot
ggplot(prot_diff_F_out, aes(x = diff, y = log_pval)) +
  geom_point(aes(color = regulation), size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  labs(
    title = paste0("Female Proteins (", upregulated_count, " Upregulated, ", downregulated_count, " Downregulated)"),
    x = "Log2 Fold Change (Difference in Abundance)",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(-2.0, 2.5) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16),  # Larger plot title
    axis.title.x = element_text(size = 14),             # Increase x-axis label size
    axis.title.y = element_text(size = 14),             # Increase y-axis label size
    axis.text = element_text(size = 12)                 # Increase axis tick text size
  ) +
  geom_text_repel(
    data = top_12_upregulated,
    aes(label = rownames(top_12_upregulated)),
    size = 6,
    box.padding = 0.4,
    max.overlaps = 100
  ) +
  geom_text_repel(
    data = top_12_downregulated,
    aes(label = rownames(top_12_downregulated)),
    size = 6,
    box.padding = 0.4,
    max.overlaps = 100
  )

# Step 1: Filter out "Negr1"
data_processed_F_out <- data_processed_F[rownames(data_processed_F) != "Negr1", ]

# Step 2: Log2 transform the data (after adding 1 to avoid -Inf)
data_pca <- t(log2(data_processed_F_out + 1))

# Step 3: Run PCA
pca_res <- prcomp(data_pca, scale. = FALSE)

# Step 4: Calculate percent variance explained
var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# Step 5: Create PCA dataframe
pca_df <- as.data.frame(pca_res$x)

# Extract group info from row names (e.g., "FKO_1", "FWT_2", etc.)
pca_df$Group <- factor(gsub("_.*", "", rownames(pca_df)))

# Step 6: Set custom colors
group_colors <- c("FKO" = "darkgreen", "FWT" = "purple")

# Step 7: Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 6, alpha = 0.8) +
  labs(
    title = "PCA of Female Protein Abundance",
    x = pc1_label,
    y = pc2_label
  ) +
  scale_color_manual(values = group_colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 12)
  )
