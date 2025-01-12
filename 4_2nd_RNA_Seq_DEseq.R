library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)

# convert transcript ID to genomic ID
tx2gene = read.table("~/Desktop/R_code/tx2gene.txt", header=TRUE)

# load sample information
samples = read.csv("~/Desktop/R_code/samples_2nd.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))

# Create a column that is 'mutant'
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))
rownames(samples) = samples$dir_name

# load sequencing reads
kallisto.base.dir = "~/Desktop/R_code/Abundance_run2"

# make empty vectors for each time point
contrast_non_shaking = c()
num_genes_non_shaking = c()
up_down_non_shaking = c()

# DEseq2 analysis
for (time in c(24, 48)){
  # CC_WT vs. AX_WT
  samples_sub = samples %>% dplyr::filter(time_hr==time & shaking!="TRUE" & cell_line == "WT")
  files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")
  txi = tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ axenic)
  
  # Filter out those genes with <10 counts in more than 1/4 of the samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  
  # set threshold for DE genes
  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_WT_CO_vs_WT_AX")); num_genes_non_shaking = append(num_genes_non_shaking, dim(up)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "up");
  
  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_WT_CO_vs_WT_AX")); num_genes_non_shaking = append(num_genes_non_shaking, dim(down)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "down");
  
  # de_genes_up_and_down_WT_non_shaking = as.data.frame(res) %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
  
  # write.csv(de_genes_up_and_down_WT_non_shaking, file.path("~/Desktop/R_code", paste0("de_genes_non_shaking_COWT_vs_AXWT", time, ".csv")))
  
  # CC_KO vs. CC_WT
  samples_sub = samples %>% dplyr::filter(time_hr==time & shaking!="TRUE" & axenic=="FALSE")
  files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")
  txi = tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ mutant)
  
  # Filter out those genes with <10 counts in more than 1/4 of the samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  
  # set threshold for DE genes
  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_mutants_vs_WT")); num_genes_non_shaking = append(num_genes_non_shaking, dim(up)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "up");
  
  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_mutants_vs_WT")); num_genes_non_shaking = append(num_genes_non_shaking, dim(down)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "down");
  
  # de_genes_up_and_down_mutant_non_shaking = as.data.frame(res) %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
  
  # write.csv(de_genes_up_and_down_mutant_non_shaking, file.path("~/Desktop/R_code", paste0("de_genes_non_shaking_mutant_vs_WT", time, ".csv")))
  
  
}

# Create the df that we will use for plotting
plotting_df_non_shaking = data.frame(contrast_non_shaking=as.factor(contrast_non_shaking), num_genes_non_shaking=num_genes_non_shaking, up_down_non_shaking=up_down_non_shaking)

# Reorder the columns so that they are in the same order as Ru's output
plotting_df_non_shaking = plotting_df_non_shaking %>% mutate(contrast_non_shaking = factor(contrast_non_shaking, levels=c("24h_WT_CO_vs_WT_AX", "24h_mutants_vs_WT", "48h_WT_CO_vs_WT_AX", "48h_mutants_vs_WT")))

# plotting
ggplot(plotting_df_non_shaking, aes(fill=up_down_non_shaking, y=num_genes_non_shaking, x=contrast_non_shaking, label=num_genes_non_shaking)) + 
  geom_bar(position="stack", stat="identity") + geom_text(size = 10, position = position_stack(vjust = 0.5)) + theme(axis.text=element_text(size=20),
                                                                                                                     axis.title=element_text(size=20), plot.title = element_text(size = 40), legend.text = element_text(size=25), axis.text.x = element_text(angle = -45, vjust=-0/5, hjust=0)) + ggtitle("non shaking") +
  scale_fill_manual(values=c("up" = "#404040", "down" = "#AFABAB"))

# ggsave("~/Desktop/R_code/DEGs.png")

# PCA plot
samples_non_shaking = samples %>% dplyr::filter(shaking!="TRUE")

# combine "time" and "co-culture" traits
samples_non_shaking = samples_non_shaking %>% mutate(treatment = as.factor(str_c(case_when(axenic == "TRUE" ~ "AX", TRUE ~ "CO"), as.character(time_hr))))

files_non_shaking <- file.path(kallisto.base.dir, samples_non_shaking$dir_name, "abundance.h5")

# read Abundunce tables and read normalization
txi_non_shaking = tximport(files_non_shaking, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_non_shaking = DESeqDataSetFromTximport(txi_non_shaking, colData = samples_non_shaking, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_non_shaking <- rowSums(counts(dds_non_shaking) >= 10) >= ceiling(dim(samples_non_shaking)[[1]]/4)
dds_non_shaking <- dds_non_shaking[keep_non_shaking,]

# Fit the model and run the DEseq2 analysis
dds_non_shaking = DESeq(dds_non_shaking)

vsd_non_shaking <- vst(dds_non_shaking, blind=FALSE)
rld_non_shaking <- rlog(dds_non_shaking, blind=FALSE)

pcaData_non_shaking = plotPCA(vsd_non_shaking, intgroup=c("cell_line", "treatment"), returnData=TRUE)

percentVar_non_shaking <- round(100 * attr(pcaData_non_shaking, "percentVar"))

ggplot(pcaData_non_shaking, aes(PC1, PC2, color=cell_line, shape=treatment)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar_non_shaking[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_non_shaking[2],"% variance")) + 
  coord_fixed() + scale_color_manual(values=c("A7" = "#082BD2", "A8" = "#65D527", "WT" = "#000000")) + ggtitle("2nd RNA-seq non-shaking PCA")+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(color = "black", fill = NA))

# ggsave("~/Desktop/R_code/PCA.png", width = 10, height = 7)

# Analyses for later WGCNA

########### AXENIC DE ###############
samples_ax_vs_co = samples %>% dplyr::filter(mutant=="FALSE" & shaking!="TRUE")

files_ax_vs_co <- file.path(kallisto.base.dir, samples_ax_vs_co$dir_name, "abundance.h5")

# read Abundunce tables and normalize
txi_ax_vs_co = tximport(files_ax_vs_co, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_ax_vs_co = DESeqDataSetFromTximport(txi_ax_vs_co, colData = samples_ax_vs_co, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_ax_vs_co <- rowSums(counts(dds_ax_vs_co) >= 10) >= ceiling(dim(samples_ax_vs_co)[[1]]/4)
dds_ax_vs_co <- dds_ax_vs_co[keep_ax_vs_co,]

# Fit the model and run the DE analysis
dds_ax_vs_co = DESeq(dds_ax_vs_co, minReplicatesForReplace=3)

res_ax_vs_co = results(dds_ax_vs_co)

res_ax_vs_co = as.data.frame(res_ax_vs_co)

# DE genes between CC_WT and AX_WT
res_ax_vs_co.filtered = res_ax_vs_co %>% dplyr::filter(padj <= 0.01) %>% dplyr::arrange(padj)
res_ax_vs_co.filtered_05 = res_ax_vs_co %>% dplyr::filter(padj <= 0.05) %>% dplyr::arrange(padj)

save(res_ax_vs_co.filtered, file="~/Desktop/R_code/run_2_res_ax_vs_co.filtered.ensembl.RData")
save(res_ax_vs_co.filtered_05, file="~/Desktop/R_code/run_2_res_ax_vs_co_05.filtered.ensembl.RData")


########### Mutant DE ###############
samples_mutant_vs_WT_co = samples %>% dplyr::filter(shaking!="TRUE" & axenic=="FALSE")

files_mutant_vs_WT_co <- file.path(kallisto.base.dir, samples_mutant_vs_WT_co$dir_name, "abundance.h5")

# read Abundunce tables and normalize
txi_mutant_vs_WT_co = tximport(files_mutant_vs_WT_co, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_mutant_vs_WT_co = DESeqDataSetFromTximport(txi_mutant_vs_WT_co, colData = samples_mutant_vs_WT_co, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_mutant_vs_WT_co <- rowSums(counts(dds_mutant_vs_WT_co) >= 10) >= ceiling(dim(samples_mutant_vs_WT_co)[[1]]/4)
dds_mutant_vs_WT_co <- dds_mutant_vs_WT_co[keep_mutant_vs_WT_co,]

# Fit the model and run the DE analysis
dds_mutant_vs_WT_co = DESeq(dds_mutant_vs_WT_co, minReplicatesForReplace=3)

res_mutant_vs_WT_co = results(dds_mutant_vs_WT_co)

res_mutant_vs_WT_co = as.data.frame(res_mutant_vs_WT_co)

# These are the DE genes for the Axenic effect
res_mutant_vs_WT_co.filtered = res_mutant_vs_WT_co %>% dplyr::filter(padj <= 0.01) %>% dplyr::arrange(padj)
res_mutant_vs_WT_co.filtered_05 = res_mutant_vs_WT_co %>% dplyr::filter(padj <= 0.05) %>% dplyr::arrange(padj)

save(res_mutant_vs_WT_co.filtered, file="~/Desktop/R_code/run_2_res_mutant_vs_WT_co.filtered.ensembl.RData")
save(res_mutant_vs_WT_co.filtered_05, file="~/Desktop/R_code/run_2_res_mutant_vs_WT_co_05.filtered.ensembl.RData")


