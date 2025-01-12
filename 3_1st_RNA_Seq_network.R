library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)
library(WGCNA)
library(gridExtra)
library(ggdendro)
library(ggrepel)
library(igraph)
library("scales")

# load in "adjacency" file generated in the WGCNA analysis
load("~/Desktop/R_code/adjacency.mat.subset.ensembl.RData")

# load de_genes
load("~/Desktop/R_code/de_genes_05.subset.ensembl.RData")

# subset the adjacency matrix to only the de_genes
adjacency_sub = adjacency[rownames(de_genes_05), rownames(de_genes_05)]

# set connection threshold
h_thresh = 0.1

adjacency_sub[adjacency_sub < h_thresh] <- 0
net <- graph_from_adjacency_matrix(adjacency_sub, mode="undirected", weighted=TRUE, diag=FALSE)
E(net)$width <- E(net)$weight * 10
V(net)$size <- 7

V(net)$color = ifelse(names(V(net)) == "Phatr3_J43365", "red", "white")

plot(net)

# Produce dataframes of the vertice and edge information 
edge_df = as_data_frame(net, what="edges")

# write.csv(edge_df, file.path("~/Desktop/R_code", paste0("DE_gene_network.csv")))

vert_df = as_data_frame(net, what="vertices")

# Load the gene significance dataframe
load("~/Desktop/R_code/geneTraitSignificance_adj.subset.ensembl.RData")

coords = layout_nicely(graph=net)
net_df_coords = data.frame(x=coords[,1], y=coords[,2], name=vert_df$name, size=vert_df$size, color=vert_df$color)
rownames(net_df_coords) = net_df_coords$name
net_df_coords$GS = geneTraitSignificance_adj[rownames(net_df_coords), "GS.axenic"]

from.x = net_df_coords$x[sapply(edge_df$from, match, net_df_coords$name)]
from.y = net_df_coords$y[sapply(edge_df$from, match, net_df_coords$name)]
to.x = net_df_coords$x[sapply(edge_df$to, match, net_df_coords$name)]
to.y = net_df_coords$y[sapply(edge_df$to, match, net_df_coords$name)]
net_edge_df = data.frame(from=edge_df$from, to=edge_df$to, from.x=from.x, from.y=from.y, to.x=to.x, to.y=to.y, weight=edge_df$weight, width=edge_df$width)

ggplot() +
  geom_segment(data=net_edge_df,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, linewidth=weight), colour="grey") +
  geom_point(data=net_df_coords, color="black", fill=net_df_coords$color, shape=21, size=rescale(net_df_coords$GS, c(10,20)), aes(x=x, y=y)) +
  geom_text(data=net_df_coords, aes(label=name, x=x, y=y)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) + xlab("") + ylab("") +
  guides(size = guide_legend("adjacency score (weight)")) +
  ggtitle("Network of DE genes (P<0.05) according to their adjacency scores Ensembl (connections >0.1 shown)")
ggsave("~/Desktop/R_code/de_adj_net.subset.ensembl.png", height=50, width=40, units="cm")