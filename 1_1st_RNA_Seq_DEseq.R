library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)
library(ggrepel)

# convert transcript ID to genomic ID
tx2gene = read.table("~/Desktop/R_code/tx2gene.txt", header=TRUE)

# load sample information
samples = read.csv("~/Desktop/R_code/samples_1st.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))

# load sequencing reads
kallisto.base.dir = "~/Desktop/R_code/Abundance_run1"

# choose time points for later DEseq2 analysis
times = c(0.5, 3, 24, 48)

# make empty vectors for collecting DEGs at each time point
contrast_time = c()
num_genes = c()
up_down = c()
de_genes = c()

# The loop for DEseq2 analysis at each time point
for (time in times){
  samples_sub = samples %>% dplyr::filter(time_hr==time)
  
  files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")
  
  # use tximport to read in the abundance tables and to normalize reads
  txi = tximport(files, type = "kallisto", tx2gene = tx2gene)
  
  # Create the DESEQ2 object
  dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ axenic)
  
  # Filter out those genes with <10 counts in more than 1/4 of the samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
  dds <- dds[keep,]
  
  # Fit the model and run the DEseq2 analysis
  dds = DESeq(dds)
  
  res = results(dds)
  
  # set threshold for DEGs
  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 1 & padj < 0.05)
  contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(up)[[1]]); up_down = append(up_down, "up");
  
  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -1 & padj < 0.05)
  contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(down)[[1]]); up_down = append(up_down, "down");
  
  # another possibility
  # de_genes_up_and_down = as.data.frame(res) %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
  # write.csv(de_genes_up_and_down, file.path("~/Desktop/R_code/R", paste0("de_genes_", time, ".csv")))
  
  # Collect the differentially expressed genes
  de_genes[[as.character(time)]] = c(rownames(up), rownames(down))
}

# plot number of DEGs
plotting_df = data.frame(contrast_time=as.factor(contrast_time), num_genes=num_genes, up_down=up_down)

ggplot(plotting_df, aes(fill=up_down, y=num_genes, x=contrast_time, label=num_genes)) + 
  geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.6)) + 
  scale_fill_manual(values=c("up" = "#404040", "down" = "#AFABAB")) + ggtitle("1st RNA-seq run DEGs")

# save figure
# ggsave("~/Desktop/R/rna1_DEGs.png")


### PCA ###
files_all_samples <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")

names(files_all_samples) = samples$dir_name

txi_all_samples = tximport(files_all_samples, type = "kallisto", tx2gene = tx2gene)

# Create the DESeq2 object
dds_all_samples = DESeqDataSetFromTximport(txi_all_samples, colData = samples, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_all_samples <- rowSums(counts(dds_all_samples) >= 10) >= ceiling(dim(samples)[[1]]/4)
dds_all_samples <- dds_all_samples[keep_all_samples,]

# Fit the model and run the DE analysis
dds_all_samples = DESeq(dds_all_samples)

vsd_all_samples <- vst(dds_all_samples, blind=FALSE)

pcaData = plotPCA(vsd_all_samples, intgroup=c("time_hr", "axenic"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time_hr, shape=axenic)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + scale_color_manual(values=c("0.5" = "#D721BB", "3" = "#144BD5", "24" = "#3CCA23", "48" = "#21CBCA")) +
  ggtitle("1st RNA-seq PCA all samples")+
  theme(text = element_text(size=10, face = "bold"),
        panel.background = element_rect(fill="transparent"),
        panel.border = element_rect(color = "black", fill = NA))

# ggsave("~/Desktop/R/rna1_PCA.png",width = 8, height = 6)

# look for common DEGs at different time points
common_genes = Reduce(intersect, de_genes)

ggvenn(
  de_genes, 
  fill_color = c("#D721BB", "#144BD5", "#3CCA23", "#21CBCA"),
  stroke_size = 0.5, set_name_size = 8
)

# ggsave("~/Desktop/R/Venn.png", width=30, height=30, unit="cm")












