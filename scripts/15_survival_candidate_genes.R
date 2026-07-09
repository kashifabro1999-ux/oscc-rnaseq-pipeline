library(dplyr)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

candidate_genes <- hub %>%
  arrange(desc(Degree)) %>%
  select(
    GeneSymbol,
    GeneName,
    Degree,
    log2FoldChange,
    padj
  )

write.table(
  candidate_genes,
  "results/survival_candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Survival candidate gene list exported successfully.\n")
