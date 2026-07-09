library(dplyr)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

biomarkers <- hub %>%
  mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated in Tumor", "Downregulated in Tumor"),
    PriorityScore = Degree * abs(log2FoldChange),
    Category = case_when(
      GeneSymbol %in% c("MMP9", "MMP14", "PLAU", "COL1A1") ~ "Invasion / ECM remodeling",
      GeneSymbol %in% c("PPARG", "ACACB", "GPD1", "LPL", "PCK1", "PLIN1", "SLC2A4") ~ "Metabolic reprogramming",
      GeneSymbol %in% c("ADH1B", "CYP3A4") ~ "Carcinogen / drug metabolism",
      TRUE ~ "Other cancer-associated hub"
    )
  ) %>%
  arrange(desc(PriorityScore))

write.table(
  biomarkers,
  "results/biomarker_priority_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Biomarker priority table created successfully.\n")
