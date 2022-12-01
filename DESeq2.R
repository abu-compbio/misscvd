# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("EnhancedVolcano")
# install.packages("tidyverse")

library(DESeq2)
library(tidyverse)

data <- read.table('RNA_count_matrix.txt')
is.na(data$x1) <- na.omit(data)
genes <- read.table('RNA_gene_metadata.txt', header = 1)
genes[genes$gene_name == 'PVRL1', 'gene_name'] = 'NECTIN1'
samples <- read.table("RNA_sample_metadata.txt", sep = ',', header=TRUE)
data <- data[!duplicated(genes$gene_name),]
genes <- genes[!duplicated(genes$gene_name),]
row.names(data) <- genes$gene_name


dds_A <- DESeqDataSetFromMatrix(countData = data[,samples[samples$condition=='A','sample_id']],
                              colData = samples[samples$condition=='A',],
                              design= ~ sex)
dds_A <- DESeq(dds_A)
resultsNames(dds_A)
res_A <- results(dds_A, contrast = c('sex', 'female', 'male'))

dds_B <- DESeqDataSetFromMatrix(countData = data[,samples[samples$condition=='B','sample_id']],
                                colData = samples[samples$condition=='B',],
                                design= ~ sex)
dds_B <- DESeq(dds_B)
resultsNames(dds_B)
res_B <- results(dds_B, contrast = c('sex', 'female', 'male'))


ernest_genes <- read.table('NGC_genes_Ernest.txt')$V1
all(ernest_genes %in% genes$gene_name) # TRUE

ernest_A <- res_A[ernest_genes,][order(res_A[ernest_genes,'padj']),]
ernest_A$padj <- p.adjust(ernest_A$pvalue, method = 'BH')
ernest_B <- res_B[ernest_genes,][order(res_B[ernest_genes,'padj']),]
ernest_B$padj <- p.adjust(ernest_B$pvalue, method = 'BH')

write.csv('DESeq2_ernest_naive.csv', x = ernest_A)
write.csv('DESeq2_ernest_IFNg.csv', x = ernest_B)


library(EnhancedVolcano)

EnhancedVolcano(ernest_A,
                lab = rownames(ernest_A),
                title = 'ernest genes naive samples (Female vs Male)',
                pCutoff = 0.05,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ italic('Padj')),
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-adj", expression(p - adj ~ and
                                                                                      ~ log[2] ~ FC))
                )

EnhancedVolcano(ernest_B,
                lab = rownames(ernest_B),
                title = 'ernest genes IFNg samples (Female vs Male)',
                pCutoff = 0.05,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ italic('Padj')),
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-adj", expression(p - adj ~ and
                                                                                    ~ log[2] ~ FC))
)




dds_male <- DESeqDataSetFromMatrix(countData = data[,samples[(samples$condition=='A' | samples$condition=='B') & samples$sex =='male','sample_id']],
                                colData = samples[(samples$condition=='A' | samples$condition=='B') & samples$sex =='male',],
                                design= ~ condition_name)
dds_male <- DESeq(dds_male)
resultsNames(dds_male)
res_male <- results(dds_male, contrast = c('condition_name', 'IFNg', 'naive'))

dds_female <- DESeqDataSetFromMatrix(countData = data[,samples[(samples$condition=='A' | samples$condition=='B') & samples$sex =='female','sample_id']],
                                   colData = samples[(samples$condition=='A' | samples$condition=='B') & samples$sex =='female',],
                                   design= ~ condition_name)
dds_female <- DESeq(dds_female)
resultsNames(dds_female)
res_female <- results(dds_female, contrast = c('condition_name', 'IFNg', 'naive'))

ernest_male   <- res_male  [ernest_genes,][order(res_male  [ernest_genes,'padj']),]
ernest_male$padj <- p.adjust(ernest_male$pvalue, method = 'BH')
ernest_female <- res_female[ernest_genes,][order(res_female[ernest_genes,'padj']),]
ernest_female$padj <- p.adjust(ernest_female$pvalue, method = 'BH')

write.csv('DESeq2_ernest_male.csv'  , x = ernest_male)
write.csv('DESeq2_ernest_female.csv', x = ernest_female)

EnhancedVolcano(ernest_male,
                lab = rownames(ernest_male),
                title = 'ernest genes male samples (IFNg vs naive)',
                pCutoff = 0.05,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ italic('Padj')),
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-adj", expression(p - adj ~ and
                                                                                    ~ log[2] ~ FC))
)

EnhancedVolcano(ernest_female,
                lab = rownames(ernest_female),
                title = 'ernest genes female samples (IFNg vs naive)',
                pCutoff = 0.05,
                FCcutoff = 1,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = bquote(~Log[2] ~ "fold change"),
                ylab = bquote(~-Log[10] ~ italic('Padj')),
                legendLabels = c("NS", expression(Log[2] ~ FC), "p-adj", expression(p - adj ~ and
                                                                                    ~ log[2] ~ FC))
)
