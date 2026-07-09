library(ggplot2)
library(dplyr)
library(tidyr)

expr <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

top_genes <- hub %>%
  arrange(desc(Degree)) %>%
  head(10)

dir.create("plots/hub_boxplots", showWarnings = FALSE)

for (i in 1:nrow(top_genes)) {
  
  entrez <- as.character(top_genes$EntrezID[i])
  symbol <- top_genes$GeneSymbol[i]
  
  if (entrez %in% rownames(expr)) {
    
    df <- data.frame(
      sample = colnames(expr),
      expression = as.numeric(expr[entrez, ])
    )
    
    df <- merge(df, sample_info, by = "sample")
    
    p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      theme_bw() +
      labs(
        title = paste0(symbol, " Expression"),
        x = "",
        y = "VST normalized expression"
      ) +
      scale_fill_manual(values = c("Control" = "#00BFC4", "Tumor" = "#F8766D")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
    
    ggsave(
      paste0("plots/hub_boxplots/", symbol, "_boxplot.png"),
      p,
      width = 5,
      height = 5,
      dpi = 300
    )
  }
}

cat("Hub gene boxplots completed successfully.\n")
