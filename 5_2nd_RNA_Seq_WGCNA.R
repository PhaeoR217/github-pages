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
library(pheatmap)
library("tidyverse")

# convert transcript ID to genomic ID
tx2gene = read.table("~/Desktop/R_code/tx2gene.txt", header=TRUE)

# load sample information
samples = read.csv("~/Desktop/R_code/samples_2nd.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))

# Create a column that is 'mutant'
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))
rownames(samples) = samples$dir_name

kallisto.base.dir = "~/Desktop/R_code/Abundance_run2"

files <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")

# read in Abundance tables and normalize
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples)[[1]]/4)
dds <- dds[keep,]

# Fit the model and run the DEseq2 analysis
dds = DESeq(dds)

# Now transform and normalize the data
# This will return log2 transformed data normalized by library size
# and with the experiment-wide mean-variance trend removed.
vsd <- vst(dds, blind=FALSE)

summary(as.data.frame(assay(vsd)))

wgcna_in = t(as.data.frame(assay(vsd)))

#==============================================
#
# Begin the WGCNA analysis
# Visualize the sample similarity
#
#==============================================
# Look for outliers
hc = hclust(dist(wgcna_in), method = "average")
ggdendrogram(hc, rotate = FALSE, size = 2)

# Check that all of the samples are good with regards to missing data
gsg = goodSamplesGenes(wgcna_in, verbose = 3)
gsg$allOK # All OK

# Create a traits file from the samples file and converted to numeric
traits = samples %>% select(axenic, time_hr, mutant) %>% mutate(axenic = as.numeric(axenic), time_hr = as.numeric(time_hr), mutant = as.numeric(mutant))

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
text_size = 5
sft_plotting_df = sft$fitIndices
r_sqr_plot = ggplot(sft_plotting_df, aes(x=Power, y=-sign(slope)*SFT.R.sq, label=Power)) + 
  geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==8, 'red', 'black')) + ggtitle("Scale independence") + 
  theme(plot.title = element_text(size = text_size, face = "bold"), 
        axis.title=element_text(size=text_size,face="bold"), 
        axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
  geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==8, 'red', 'black')) + ggtitle("Mean connectivity") + 
  theme(plot.title = element_text(size = text_size, face = "bold"),
        axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
  xlab("Soft Threshold (Power)")

threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)

# ggsave("~/Desktop/R_code/run_2_wgcna_threshold_selection.ensembl.png", plot=threshold_plot, width=20, height=10, units="cm")

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40)

if(file.exists("~/Desktop/R_code/adjacency_TOM.ensembl.RData")){
  load("~/Desktop/R_code/adjacency_TOM.ensembl.RData")
}else{
  # Make adjacency matrix
  adjacency = adjacency(wgcna_in, power = 8)
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  save(adjacency, TOM, file="~/Desktop/R_code/adjacency_TOM.ensembl.RData")
}

dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
ggdendrogram(geneTree, rotate = FALSE, size = 2)

# set the minimum module size
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# 0 is automatically assigned to grey in the labels2colors method.
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

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

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(wgcna_in, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(200));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#==============================================
#
# Module-trait analysis
#
#==============================================
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

png(filename="~/Desktop/R_code/mod_trait.ensembl.png", width=10, height=10, units="cm", res=300)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
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
# Define variable "axenic"
axenic = as.data.frame(traits$axenic)
names(axenic) = "axenic"
geneTraitSignificance_axenic = as.data.frame(cor(wgcna_in, axenic, use = "p"))
GSPvalue_axenic = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_axenic), nSamples))
names(geneTraitSignificance_axenic) = paste("GS.", names(axenic), sep="")
names(GSPvalue_axenic) = paste("p.GS.", names(axenic), sep="")

# Define variable "time"
time_hr = as.data.frame(traits$time_hr)
names(time_hr) = "time_hr"
geneTraitSignificance_time_hr = as.data.frame(cor(wgcna_in, time_hr, use = "p"))
GSPvalue_time_hr = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_time_hr), nSamples))
names(geneTraitSignificance_time_hr) = paste("GS.", names(time_hr), sep="")
names(GSPvalue_time_hr) = paste("p.GS.", names(time_hr), sep="")

# Define variable "mutant"
mutant = as.data.frame(traits$mutant)
names(mutant) = "mutant"
geneTraitSignificance_mutant = as.data.frame(cor(wgcna_in, mutant, use = "p"))
GSPvalue_mutant = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mutant), nSamples))
names(geneTraitSignificance_mutant) = paste("GS.", names(mutant), sep="")
names(GSPvalue_mutant) = paste("p.GS.", names(mutant), sep="")

# Put all four of the dataframes together into a single df
l = list(geneTraitSignificance_axenic, GSPvalue_axenic, geneTraitSignificance_time_hr, GSPvalue_time_hr, geneTraitSignificance_mutant, GSPvalue_mutant)
geneTraitSignificance = Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(geneTraitSignificance) = geneTraitSignificance$rn
geneTraitSignificance = geneTraitSignificance %>% select(-rn)

# write.csv(geneTraitSignificance_mutant, file.path("~/Desktop/R_code", paste0("Gene_trait_significance_mutant_all_samples.csv")))
# write.csv(GSPvalue_mutant, file.path("~/Desktop/R_code", paste0("Gene_trait_significance_mutant_pvalue_all_samples.csv")))

#==============================================
#
# Intra module investigation
#
#==============================================
# select a module of interest
module_of_interest = "cyan"
  column = match(module_of_interest, modNames)
  moduleGenes = moduleColors==module_of_interest
  
  plotting.df = data.frame(MM=abs(geneModuleMembership[moduleGenes, column]), GS=abs(geneTraitSignificance[moduleGenes, "GS.mutant"]))
  row.names(plotting.df) = row.names(geneModuleMembership)[moduleGenes]
  plotting.df$gene_names = row.names(geneModuleMembership)[moduleGenes]
  
  ggplot(plotting.df, aes(x=MM, y=GS, label=gene_names)) + 
    geom_point(color=dplyr::case_when(plotting.df$GS > 0.9 ~ "red", TRUE ~ "black")) + 
    geom_label_repel(aes(label=ifelse(GS > 0.9, gene_names, ""))) + 
    xlab(paste("Module Membership in", module_of_interest, "module")) + ylab("Gene significance for axenic") + 
    ggtitle("Module membership vs. gene significance\n")
  #ggsave("~/Desktop/R_code/cyan_MM_GS_scatter.ensembl.png")

geneModuleMembership[head(rownames(geneTraitSignificance %>% arrange(desc(abs(GS.mutant)))), 100),]

gs_mut = (geneTraitSignificance %>% dplyr::filter(GS.mutant>0.5))$GS.mutant
hist_df = data.frame(x=gs_mut)
ggplot(hist_df, aes(x=x)) + geom_histogram()  

module_gene_names = rownames(geneModuleMembership)[moduleGenes]

save(geneTraitSignificance, file="~/Desktop/R_code/geneTraitSignificance.non_shaking.ensembl.RData")
save(module_gene_names, file="~/Desktop/R_code/cyan_module_genes_names.non_shaking.ensembl.RData")
save(adjacency, file="~/Desktop/R_code/adjacency.non_shaking.ensembl.RData")
# save.image("~/Desktop/R_code/run_2_wgcna.ensembl.RData")

# "mutant" related network
load("~/Desktop/R_code/run_2_res_mutant_vs_WT_co_05.filtered.ensembl.RData")
res_mutant_vs_WT_co.filtered = res_mutant_vs_WT_co.filtered_05

head(res_mutant_vs_WT_co.filtered)
dim(res_mutant_vs_WT_co.filtered)
sum(module_gene_names %in% rownames(res_mutant_vs_WT_co.filtered))
length(module_gene_names)

res_mutant_vs_WT_co.filtered$loc = 1
res_mutant_vs_WT_co.filtered$rn = rownames(res_mutant_vs_WT_co.filtered)
ggplot(res_mutant_vs_WT_co.filtered, aes(x=loc, y=-log10(padj))) + 
  geom_point(color=dplyr::case_when(res_mutant_vs_WT_co.filtered$rn %in% module_gene_names ~ "red", TRUE ~ "black"), position = "jitter") + 
  ggtitle("Mutant sig genes (red== in cyan module membership)")

# plot up the relationship between GS.mutant and padj of the DE genes
dim(res_mutant_vs_WT_co.filtered)
res_mutant_vs_WT_co.filtered.nomissing = res_mutant_vs_WT_co.filtered[rownames(res_mutant_vs_WT_co.filtered) %in% rownames(geneTraitSignificance),]
dim(res_mutant_vs_WT_co.filtered.nomissing)

res_mutant_vs_WT_co.filtered.nomissing = res_mutant_vs_WT_co.filtered[rownames(res_mutant_vs_WT_co.filtered) %in% rownames(geneTraitSignificance),]

res_mutant_vs_WT_co.filtered.nomissing$GS.mutant = geneTraitSignificance[rownames(res_mutant_vs_WT_co.filtered.nomissing), "GS.mutant"]

ggplot(res_mutant_vs_WT_co.filtered.nomissing, aes(x=-log10(padj), y=abs(GS.mutant))) + geom_point() +
  geom_point(color=dplyr::case_when(res_mutant_vs_WT_co.filtered.nomissing$rn %in% module_gene_names ~ "red", TRUE ~ "black"), position = "jitter") +
  scale_x_continuous(trans='log10') + ggtitle("padj vs GS.mutant for the mutant DE genes (red=cyan module membership)")

de_gene_names = rownames(res_mutant_vs_WT_co.filtered.nomissing)
adjacency_de = adjacency[de_gene_names, de_gene_names]

save(adjacency_de, file="~/Desktop/R_code/adjacency_de.non_shaking.ensembl.RData")

adjacency_df = data.frame(adj_sum=rowSums(adjacency_de), lc_member=rownames(adjacency_de) %in% module_gene_names)
ggplot(adjacency_df, aes(x=lc_member, y=adj_sum)) + geom_violin(trim=FALSE) + geom_point() + ggtitle("adjacency_score_vs_tan_membership")

load("~/Desktop/R_code/cyan_module_genes_names.non_shaking.ensembl.RData")

# To identify the other modules we should create a tree and then perform hierarchical clustering on that tree
hc <- hclust(dist(1-adjacency_de), method = "complete")
as.dendrogram(hc) %>% plot(horiz = TRUE)
abline(v=5.8, col = "red")
gene_clusters = cutree(hc, k = 20)
gene_clusters_df = as.data.frame(gene_clusters, row.names=names(gene_clusters))

dev.off()

# To look up the genes in a particular cluster (cluster==12):
# cluster = 12
# Gene_names_cluster = rownames(gene_clusters_df %>% dplyr::filter(gene_clusters==cluster))
# Gene_names_cluster_num = as.data.frame(Gene_names_cluster)
# write.csv(Gene_names_cluster_num, file.path("~/Desktop/R_code", paste0("Cluster_", cluster, "_gene_names.csv")))

cyan_membership_df = as.data.frame(ifelse(rownames(adjacency_de) %in% module_gene_names, "yes", "no"), row.names=rownames(adjacency_de), col.names=c("tan_m"))
colnames(cyan_membership_df) = "tan_m"
row_annotations = merge(gene_clusters_df, cyan_membership_df, by="row.names")
rownames(row_annotations) = row_annotations$Row.names

row_annotations = row_annotations %>% mutate(gene_clusters=as.factor(gene_clusters)) %>% 
  mutate(cluster_1=as.factor(gene_clusters == 1)) %>% 
  mutate(cluster_9=as.factor(gene_clusters == 3)) %>% 
  mutate(cluster_11_13=as.factor(gene_clusters == 9 | gene_clusters == 12 | gene_clusters == 7)) %>% 
  dplyr::select(-Row.names)
png("~/Desktop/R_code/mutant.de.non_shaking.adjacency.dendro.potential_networks.ensembl.png", width=40, height=30, units="cm", res=300)

pheatmap(as.dist(1-adjacency_de), annotation_row=row_annotations, main="DE genes for mutant trait, clustered by adjacency matrix")

dev.off()

# save(gene_clusters, file="~/Desktop/R_code/mutant.de.non_shaking.adjacency.dendro.clusters.ensembl.RData")


# "axenic" related network
load("~/Desktop/R_code/run_2_res_ax_vs_co_05.filtered.ensembl.RData")
load("~/Desktop/R_code/geneTraitSignificance.non_shaking.ensembl.RData")
res_ax_vs_co.filtered = res_ax_vs_co.filtered_05

head(res_ax_vs_co.filtered)
head(geneTraitSignificance)
dim(res_ax_vs_co.filtered)

res_ax_vs_co.filtered.nomissing = res_ax_vs_co.filtered[rownames(res_ax_vs_co.filtered) %in% rownames(geneTraitSignificance),]
dim(res_ax_vs_co.filtered.nomissing)

de_gene_names_axenic = rownames(res_ax_vs_co.filtered.nomissing)
adjacency_de_axenic = adjacency[de_gene_names_axenic, de_gene_names_axenic]

save(adjacency_de_axenic, file="~/Desktop/R_code/adjacency_de_axenic.non_shaking.ensembl.RData")

annot_df = data.frame(padj=-log10(res_ax_vs_co.filtered[rownames(adjacency_de_axenic),]$padj), abs_lg2fc=abs(res_ax_vs_co.filtered[rownames(adjacency_de_axenic),]$log2FoldChange), row.names=rownames(adjacency_de_axenic))

png("~/Desktop/R_code/temp_plot.ensembl.png", width=40, height=30, units="cm", res=300)

pheatmap(as.dist(1-adjacency_de_axenic), annotation_row=annot_df)
dev.off()

# identify the other modules
hc_axenic <- hclust(dist(1-adjacency_de_axenic), method = "complete")
png("~/Desktop/R_code/temp_plot.ensembl.png", width=40, height=30, units="cm", res=300)

as.dendrogram(hc_axenic) %>% plot(horiz = TRUE)
abline(v=4.8, col = "red")
dev.off()

gene_clusters_axenic = cutree(tree = hc_axenic, k = 20)
gene_clusters_df_axenic = data.frame(gene_clusters = gene_clusters_axenic, row.names=names(gene_clusters_axenic))
gene_clusters_df_axenic = gene_clusters_df_axenic %>% mutate(gene_clusters=as.factor(gene_clusters))

# > table(gene_clusters_axenic)

png("~/Desktop/R_code/axenic.de.non_shaking.adjacency.dendro.potential_networks.ensembl.png", width=40, height=30, units="cm", res=300)

pheatmap(as.dist(1-adjacency_de_axenic), annotation_row=annot_df, annotation_col=gene_clusters_df_axenic, main="DE genes for axenic trait, clustered by adjacency matrix")

dev.off()

# save.image("~/Desktop/R_code/run_2_wgcna.ensembl.RData")
