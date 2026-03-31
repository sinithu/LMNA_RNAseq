#Differential gene expression analysis

#Reading the count matrix
count_matrix <- read.delim("count_matrix.txt")
count_matrix <- as.data.frame(count_matrix) #columns: gene_id, C1WT, C2WT, C3WT, LMN1, LMN2, LMN3, gene_name, gene_chr, gene_start, gene_end, gene_strand, gene_length, gene_biotype, gene_description, tf_family

#Modifying the matrix
rownames(count_matrix) <- count_matrix$gene_id #gene names as rownames
count_matrix_protein_coding <- count_matrix[count_matrix$gene_biotype=="protein_coding",] #including only protein coding genes
cts <- count_matrix_protein_coding[,2:7] #including only sample read count columns

#Annotating samples
coldata <- as.data.frame(matrix(c(colnames(cts), rep('healthy', 3), rep('disease', 3)), nrow = 6, ncol = 2))
colnames(coldata) <- c('id', 'condition')
coldata$condition <- factor(coldata$condition, levels = c("healthy", "disease")) #defining the reference level

#DGE analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata,design= ~ condition)
dds <- DESeq(dds)

#Results
res05 <- results(dds, alpha=0.05)
res05_ordered <- res05[order(res05$padj),] #ordering the results by adjusted p-value
filtered_res05 <- subset(res05_ordered, padj < 0.05 & abs(log2FoldChange) > 1)
filtered_res05_up <- subset(res05_ordered, padj < 0.05 & log2FoldChange > 1)
filtered_res05_down <- subset(res05_ordered, padj < 0.05 & log2FoldChange < -1)