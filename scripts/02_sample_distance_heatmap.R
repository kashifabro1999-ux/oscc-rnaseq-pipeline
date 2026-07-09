library(pheatmap)
library(RColorBrewer)

vst <- read.delim("data/vst_counts.tsv",
                  row.names = 1,
                  check.names = FALSE)

sample_info <- read.delim("data/sample_info.tsv",
                          check.names = FALSE)

vst <- as.matrix(vst)

sample_dists <- dist(t(vst))

sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- colnames(vst)
colnames(sample_dist_matrix) <- colnames(vst)

annotation <- data.frame(condition = sample_info$condition)
rownames(annotation) <- sample_info$sample

pheatmap(sample_dist_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(
           rev(brewer.pal(9, "Blues"))
         )(255),
         main = "Sample-to-Sample Distance Heatmap",
         filename = "plots/sample_distance_heatmap.png",
         width = 10,
         height = 10)

pdf("plots/sample_distance_heatmap.pdf",
    width = 10,
    height = 10)

pheatmap(sample_dist_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(
           rev(brewer.pal(9, "Blues"))
         )(255),
         main = "Sample-to-Sample Distance Heatmap")

dev.off()
