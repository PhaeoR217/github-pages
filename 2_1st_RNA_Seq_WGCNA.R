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
library(tidyr)

# choose samples
sample_subset = c(0.5, 3, 24, 48)

# convert transcript ID to genomic ID
tx2gene = read.table("~/Desktop/R_code/tx2gene.txt", header=TRUE)

# load sample information
samples = read.csv("~/Desktop/R_code/samples_1st.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))

# subset the samples to the relevant subset
samples = samples %>% dplyr::filter(time_hr %in% sample_subset)

rownames(samples) = lapply(lapply(str_split(samples$dir_name, "_"), FUN=head, -4), str_c, collapse="_")

# load sequencing reads
kallisto.base.dir = "~/Desktop/R_code/Abundance_run1"

files <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples)[[1]]/4)
dds <- dds[keep,]

# Fit the model and run the DE analysis
dds = DESeq(dds)

vsd <- vst(dds, blind=FALSE)

wgcna_in = t(as.data.frame(assay(vsd)))

# output csv of vst transformed and normalized samples
# write.csv(as.data.frame(assay(vsd)), "~/Desktop/R_code/vsd.counts.csv")

#==============================================
#
# Begin the WGCNA analysis
# Visualize the sample similarity
#
#==============================================
# Look for outliers
hc = hclust(dist(wgcna_in), method = "average")
#ggdendrogram(hc, rotate = FALSE, size = 2)
dd_figure <- ggdendrogram(hc, rotate = FALSE, size = 2)
dd_figure$layers[[2]]$aes_params$size <- 1
dd_figure$layers[[2]]$aes_params$colour <- "black"
  text_format <- element_text(face = "bold", color = "black", size = 15)
  dd_figure + theme(axis.text = text_format) + scale_x_discrete(limit = c('CC_48h_1',	'CC_48h_2',	'CC_48h_3','AX_48h_3','AX_48h_1','AX_48h_2',
                                                                          'CC_3h_3','CC_3h_1',	'CC_3h_2','AX_3h_1',	'AX_3h_2',	'AX_3h_3',
                                                                          'CC_24h_1',	'CC_24h_2',	'CC_24h_3','AX_24h_3','AX_24h_1',	'AX_24h_2',	'CC_0.5h_2',
                                                                          'CC_0.5h_3',	'AX_0.5h_3',	'AX_0.5h_1',	'AX_0.5h_2'
  ))

# Check that all of the samples are good with regards to missing data
gsg = goodSamplesGenes(wgcna_in, verbose = 3)
gsg$allOK

# Create a traits file from the samples file that contains the factors we are interested in (axenic, time_hr)
# converted to numeric
traits = samples %>% select(axenic, time_hr) %>% mutate(axenic = as.numeric(axenic), 
                                                        time_hr = as.numeric(time_hr))

traitColors = numbers2colors(traits, signed = FALSE)

plotDendroAndColors(hc, traitColors,
                    groupLabels = names(traits), 
                    main = "Sample dendrogram and trait heatmap")

#==============================================
#
# Construct WGCNA network
#
#==============================================

###############################################
# CHOOSING SOFT THRESHOLD POWER
###############################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_in, powerVector = powers, verbose = 5, blockSize=20000)

# Plot scale-free toplogy model fit R^2 value against power
text_size = 15
sft_plotting_df = sft$fitIndices
r_sqr_plot = ggplot(sft_plotting_df, aes(x=Power, y=-sign(slope)*SFT.R.sq, label=Power)) + 
  geom_point(size=1) + geom_text(hjust=0, vjust=0, size=5, color=ifelse(sft_plotting_df$Power==12, 'red', 'black')) + ggtitle("Scale independence") + 
  theme(plot.title = element_text(size = text_size, face = "bold"), 
        axis.title=element_text(size=text_size,face="bold"), 
        axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
  geom_point(size=1) + geom_text(hjust=0, vjust=0, size=5, color=ifelse(sft_plotting_df$Power==12, 'red', 'black')) + ggtitle("Mean connectivity") + 
  theme(plot.title = element_text(size = text_size, face = "bold"),
        axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
  xlab("Soft Threshold (Power)")

threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)

# ggsave("~/Desktop/R_code/run_1_wgcna_threshold_selection_all.png", plot=threshold_plot, width=20, height=10, units="cm")

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40)

if(file.exists("~/Desktop/R_code/adjacency_TOM_3_24_48h.RData")){
  load("~/Desktop/R_code/adjacency_TOM_3_24_48h.RData")
}else{
  # Make adjacency matrix
  adjacency = adjacency(wgcna_in, power = 12)
  # Turn adjacency into topological overlap. This step takes time
  TOM = TOMsimilarity(adjacency)
  save(adjacency, TOM, file="~/Desktop/R_code/adjacency_TOM_3_24_48h.RData")
}

dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# ggdendrogram(geneTree, rotate = FALSE, size = 2)

# Set the minimum module size
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# 0 is automatically assigned to grey by default
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Run the module-trait analysis on the non-merged modules
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, colors = dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

png("~/Desktop/R_code/run_1_module_trait.unmerged.subset.ensembl.png", width=20, height=60, units="cm", res=300)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships (unmerged)"))

dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(wgcna_in, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25

# Plot the cut line
abline(h=MEDissThres, col = "red")

# Merge modules
merge = mergeCloseModules(wgcna_in, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;
table(mergedColors)

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
mergedMEDiss = 1-cor(mergedMEs)
mergedMETree = hclust(as.dist(mergedMEDiss), method = "average")
plot(mergedMETree, main = "Clustering of merged module eigengenes",
     xlab = "", sub = "", cex = 1, cex.main = 1, cex.lab = 1, cex.axis = 1)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(200));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
mergedMEDiss = 1-cor(mergedMEs)
mergedMETree = hclust(as.dist(mergedMEDiss), method = "average")
plot(mergedMETree, main = "Clustering of merged module eigengenes",
     xlab = "", sub = "", cex = 1.5, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)

gene_module = as.data.frame(moduleColors)

dev.off()

# write.csv(gene_module, file.path("~/Desktop/R_code", paste0("Gene_module_info.csv")))

#==============================================
#
# Module-trait analysis
#
#==============================================

# For the time being we will continue with all of the samples
# Define numbers of genes and samples
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

png("~/Desktop/R_code/run_1_module_trait.merged.subset.ensembl.png", width=20, height=30, units="cm", res=300)

# fit the heatmap to screen size
par(mar = c(7,12,3,1))

# generating the heatmap
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =  1.2,
               cex.lab.x = 1.5,
               xLabelsAngle = 0,
               cex.lab.y = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships (merged)"))
dev.off()

#==============================================
#
# Gene Module Membership
#
#==============================================

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(wgcna_in, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#==============================================
#
# Gene Trait Significance
#
#==============================================

# Define variable axenic
axenic = as.data.frame(traits$axenic)
names(axenic) = "axenic"
geneTraitSignificance_axenic = as.data.frame(cor(wgcna_in, axenic, use = "p"))
GSPvalue_axenic = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_axenic), nSamples))
names(geneTraitSignificance_axenic) = paste("GS.", names(axenic), sep="")
names(GSPvalue_axenic) = paste("p.GS.", names(axenic), sep="")

# Define variable time_hr
time_hr = as.data.frame(traits$time_hr)
names(time_hr) = "time_hr"
geneTraitSignificance_time_hr = as.data.frame(cor(wgcna_in, time_hr, use = "p"))
GSPvalue_time_hr = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_time_hr), nSamples))
names(geneTraitSignificance_time_hr) = paste("GS.", names(time_hr), sep="")
names(GSPvalue_time_hr) = paste("p.GS.", names(time_hr), sep="")

# Put all four of the dataframes together into a single dataframe
l = list(geneTraitSignificance_axenic, GSPvalue_axenic, geneTraitSignificance_time_hr, GSPvalue_time_hr)
geneTraitSignificance = Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(geneTraitSignificance) = geneTraitSignificance$rn
geneTraitSignificance = geneTraitSignificance %>% select(-rn)

# write.csv(geneTraitSignificance, file.path("~/Desktop/R_code", paste0("Gene_trait_significance_run1_3_24_48h.csv")))
# write.csv(geneModuleMembership, file.path("~/Desktop/R_code", paste0("Gene_module_membership_run1_3_24_48h.csv")))
# write.csv(MMPvalue, file.path("~/Desktop/R_code", paste0("Gene_module_membership_pvalue_run1_3_24_48h.csv")))

# Search for the top 10 genes that relate to axenic state
head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)

# Get a connectivity score for each of the genes
adjacency_sum = rowSums(adjacency)
temp_merged = merge(geneTraitSignificance, adjacency_sum, by="row.names")
rownames(temp_merged) = temp_merged$Row.names
geneTraitSignificance_adj = temp_merged %>% mutate(adj_score=y) %>% select(-Row.names, -y) %>% arrange(desc(abs(GS.axenic)))
head(geneTraitSignificance_adj)
head(geneTraitSignificance_adj,10)

# Phatr3_J43365 has high gene significance
gene_43365 = "Phatr3_J43365"
geneTraitSignificance_adj[gene_43365,]

geneModuleMembership[gene_43365,]


# Look for genes highly connected to Phatr3_J43365
adj_plot_df = data.frame(adj=as.numeric(adjacency[gene_43365,]))
row.names(adj_plot_df) = names(adjacency[gene_43365,])
summary(adj_plot_df$adj)
head(adj_plot_df %>% dplyr::arrange(desc(adj)), 20)
adj_plot_df = adj_plot_df %>% mutate(gene = gene_43365, label=row.names(adj_plot_df))
# ggplot(adj_plot_df, aes(x=gene, y=adj, lab=label)) + geom_point() + geom_label_repel(aes(label=ifelse(adj > .1, label, ""))) + ggtitle(paste0("Adjacency score to ", gene_43365))

save(adj_plot_df, file="~/Desktop/R_code/adj_plot_df.ensembl.RData")

get_de_network_size = function(gene_in_q){
  adj_plot_df_gene_in_q = data.frame(adj=as.numeric(adjacency[gene_in_q,]))
  row.names(adj_plot_df_gene_in_q) = names(adjacency[gene_in_q,])
  adj_plot_df_gene_in_q = adj_plot_df_gene_in_q %>% mutate(gene = gene_in_q, label=row.names(adj_plot_df_gene_in_q)) %>% arrange(desc(adj)) %>% dplyr::filter(adj>0.1)
  
  list(sum(row.names(adj_plot_df_gene_in_q) %in% de_genes_names), row.names(adj_plot_df_gene_in_q)[row.names(adj_plot_df_gene_in_q) %in% de_genes_names])
}

res = results(dds)
res.df = as.data.frame(res) %>% dplyr::arrange(padj)
de_genes = res.df %>% dplyr::filter(padj<0.01) %>% arrange(padj)
de_genes_05 = res.df %>% dplyr::filter(padj<0.05) %>% arrange(padj)
de_genes_names = row.names(de_genes)

net_scores = lapply(as.list(de_genes_names), FUN=get_de_network_size)
names(net_scores) = de_genes_names

save(adjacency, file="~/Desktop/R_code/adjacency.mat.subset.ensembl.RData")

save(geneTraitSignificance_adj, file="~/Desktop/R_code/geneTraitSignificance_adj.subset.ensembl.RData")

save(de_genes_05, file="~/Desktop/R_code/de_genes_05.subset.ensembl.RData")

# check gene and module relationship
names(moduleColors) = colnames(wgcna_in)

head(moduleColors)

geneTraitSignificance_adj_col = (merge(geneTraitSignificance_adj, moduleColors, by="row.names"))
rownames(geneTraitSignificance_adj_col) = geneTraitSignificance_adj_col$Row.names

geneTraitSignificance_adj_col = geneTraitSignificance_adj_col %>% arrange(desc(abs(GS.axenic))) %>% select(-Row.names)

table((geneTraitSignificance_adj_col %>% dplyr::filter(abs(GS.axenic) > 0.7) %>% dplyr::select(y))$y)
head(geneTraitSignificance_adj_col, 40)

head(res.df)
genes_for_plot = rownames(res.df)[rownames(res.df) %in% rownames(geneTraitSignificance_adj_col)]

de_padj_plot_df = data.frame(GS.axenic=geneTraitSignificance[genes_for_plot, "GS.axenic"], padj=res.df[genes_for_plot, "padj"])
rownames(de_padj_plot_df) = genes_for_plot
ggplot(de_padj_plot_df, aes(x=abs(GS.axenic), y=-log10(padj))) + geom_point(color=ifelse(de_padj_plot_df$padj < 0.05, "red", "black"))

# ggsave("~/Desktop/R_code/rna_1.padj.GS.axenic.subset.ensembl.png")

# correlate gene significance to module membership
geneModuleMembership[gene_43365,]
# gene_43365 is in module grey

module_of_interest = "grey"
  column = match(module_of_interest, modNames)
  moduleGenes = moduleColors==module_of_interest
  moduleGenes_names = row.names(geneModuleMembership)[moduleGenes]
  plotting.df = data.frame(MM=abs(geneModuleMembership[moduleGenes, column]), GS=abs(geneTraitSignificance[moduleGenes, "GS.axenic"]), padj=res.df[moduleGenes_names, "padj"])
  row.names(plotting.df) = row.names(geneModuleMembership)[moduleGenes]
  plotting.df$gene_names = row.names(geneModuleMembership)[moduleGenes]

  # color the genes according to their Padj value
  ggplot(plotting.df, aes(x=MM, y=GS, label=gene_names)) + 
    geom_point(color=dplyr::case_when(plotting.df$padj < 0.05 ~ "red", TRUE ~ "black")) + 
    geom_label_repel(aes(label=ifelse(GS > .75, gene_names, ""))) + 
    xlab(paste("Module Membership in", module_of_interest, "module")) + ylab("Gene correlation for axenic state") + 
    ggtitle("Module membership vs. gene significance; red= DE padj<0.05\n")
  
# ggsave("~/Desktop/R_code/skyblue1_MM_GS_scatter.subset.ensembl.png")
  
# save.image("~/Desktop/R_code/wgcna.subset.ensembl.RData")