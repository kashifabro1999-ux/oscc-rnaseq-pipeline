library(igraph)
library(ggraph)
library(ggplot2)

nodes <- read.delim("results/PPI_nodes_top300_DEGs.tsv", check.names = FALSE)
edges <- read.delim("results/PPI_edges_top300_DEGs.tsv", check.names = FALSE)

network <- graph_from_data_frame(
  d = edges[, c("from", "to")],
  vertices = nodes,
  directed = FALSE
)

V(network)$degree <- degree(network)
V(network)$logFC <- as.numeric(V(network)$log2FoldChange)
V(network)$gene <- V(network)$GeneSymbol

png("plots/STRING_network_top300.png",
    width = 2600,
    height = 2400,
    res = 300)

ggraph(network, layout = "fr") +
  geom_edge_link(alpha = 0.25, colour = "grey70") +
  geom_node_point(aes(size = degree, color = logFC)) +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(name = "Degree") +
  theme_void() +
  ggtitle("Top 300 DEG Protein-Protein Interaction Network")

dev.off()

pdf("plots/STRING_network_top300.pdf",
    width = 12,
    height = 10)

ggraph(network, layout = "fr") +
  geom_edge_link(alpha = 0.25, colour = "grey70") +
  geom_node_point(aes(size = degree, color = logFC)) +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(name = "Degree") +
  theme_void() +
  ggtitle("Top 300 DEG Protein-Protein Interaction Network")

dev.off()

cat("STRING network visualization completed successfully.\n")
