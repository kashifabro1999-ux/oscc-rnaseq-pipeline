library(ggplot2)

vst <- read.delim("data/vst_counts.tsv", row.names = 1, check.names = FALSE)
sample_info <- read.delim("data/sample_info.tsv", check.names = FALSE)

vst <- as.matrix(vst)

pca <- prcomp(t(vst), scale. = FALSE)

pca_data <- as.data.frame(pca$x)
pca_data$sample <- rownames(pca_data)

pca_data <- merge(pca_data, sample_info, by = "sample")

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 2.5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot: Tumor vs Control")

ggsave("plots/PCA_Tumor_vs_Control.png", p, width = 8, height = 6, dpi = 300)
ggsave("plots/PCA_Tumor_vs_Control.pdf", p, width = 8, height = 6)

print(p)
