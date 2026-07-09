library(STRINGdb)
library(dplyr)

options(timeout = 600)

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

res <- res[!is.na(res$padj), ]
res <- res[!is.na(res$GeneSymbol), ]

sig <- res[
  res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]

top_genes <- sig %>%
  arrange(padj) %>%
  head(300)

dir.create("stringdb_cache", showWarnings = FALSE)

string_db <- STRINGdb$new(
  version = "12",
  species = 9606,
  score_threshold = 400,
  input_directory = "stringdb_cache"
)

mapped <- string_db$map(
  top_genes,
  "GeneSymbol",
  removeUnmappedRows = TRUE
)

interactions <- string_db$get_interactions(mapped$STRING_id)

interactions_filtered <- interactions[
  interactions$from %in% mapped$STRING_id &
    interactions$to %in% mapped$STRING_id,
]

nodes <- mapped[, c(
  "STRING_id",
  "EntrezID",
  "GeneSymbol",
  "GeneName",
  "log2FoldChange",
  "padj"
)]

degree_table <- data.frame(
  STRING_id = c(interactions_filtered$from, interactions_filtered$to)
) %>%
  group_by(STRING_id) %>%
  summarise(Degree = n(), .groups = "drop") %>%
  arrange(desc(Degree))

hub_genes <- merge(
  degree_table,
  nodes,
  by = "STRING_id"
)

hub_genes <- hub_genes[order(-hub_genes$Degree), ]

write.table(
  nodes,
  "results/PPI_nodes_top300_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interactions_filtered,
  "results/PPI_edges_top300_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hub_genes,
  "results/hub_genes_ranked.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  head(hub_genes, 30),
  "results/top30_hub_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("PPI and hub gene analysis completed successfully.\n")
cat("Significant DEGs:", nrow(sig), "\n")
cat("Top genes used:", nrow(top_genes), "\n")
cat("Mapped genes:", nrow(mapped), "\n")
cat("Interactions:", nrow(interactions_filtered), "\n")
cat("Top hub genes saved in results/top30_hub_genes.tsv\n")
