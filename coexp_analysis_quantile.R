library(pheatmap)
library(ggplot2)
library(WGCNA)
library(DESeq2)
library("RColorBrewer")
library(openxlsx)

cluster1 <- c('BOC',  'EPHB2',  'NPTXR',  'NDNF',  'NFASC',  'NTN5',  'SEMA3B',  'DSCAM',  'SEMA3A',  'ROBO1',  'DCC',  'EPHA5',  'RELN',  'SEMA3D',  'GFRA1',  'IGSF9',  'PRTG',  'ROBO2',  'SEMA6D',  'RGMA',  'NTRK3',  'FLRT2',  'FZD3',  'BMP7',  'NELL1',  'EPHA10',  'SEMA3E',  'PTK7',  'SEMA4C',  'SEMA4F',  'DRAXIN',  'EPHA1',  'EPHA3',  'CNTN4',  'NPFFR2',  'SEMA6A',  'BDNF',  'EPHA4',  'NTRK2',  'NRG1',  'WNT5A',  'EPHA6',  'CNTN6',  'EPHA7',  'NTNG2',  'PLXNA4',  'UNC5C',  'FGF20',  'FLRT3',  'NRCAM',  'SLIT1',  'EPHB1',  'NTRK1',  'CDH4',  'CDNF',  'SEMA5A',  'UNC5D',  'PDGFRA',  'WNT3',  'GDNF',  'PLXNA3',  'NRP1',  'SHH',  'EPHB4',  'EFNA5',  'LRTM1',  'NTNG1')
cluster2 <- c('GDF15',  'SEMA3C',  'SEMA4B',  'NOTCH2',  'SEMA4D',  'EPHA2',  'UNC5B',  'EDN1',  'EPHB6',  'SEMA3F',  'PLXNB1',  'NOTCH3',  'PDGFRB',  'ALCAM',  'NPTN',  'NENF',  'CXCL12',  'PLXNB2',  'ARRB2',  'PLXNA1',  'CSF1R',  'PLXND1',  'NOTCH1',  'ROBO3',  'BSG',  'LRP1',  'MANF')
cluster3 <- c('FRS2',  'EFNB2',  'SLIT2',  'GDF7',  'NTN4',  'PLXNC1',  'SEMA4A',  'NECTIN1',  'PTPRO',  'ROBO4',  'SEMA5B',  'NGFR',  'PLXNA2',  'NTF3',  'NTN1',  'EFNB1',  'NOG',  'PLXNB3',  'EPHB3',  'NEO1',  'ADAMTSL1',  'SEMA3G',  'LGR6',  'SLIT3',  'ARRB1',  'NRP2',  'GFRA2',  'SCN1B',  'SEMA6B')



dds_depth <- DESeqDataSetFromMatrix(counts, colData=data.frame(sample=colnames(counts)), ~1)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)


data <- read.table('./input_data/RNA_cqn_matrix.txt')
counts = data
genes <- read.table('./input_data/RNA_gene_metadata.txt', header = 1)
samples <- read.table("./input_data/RNA_sample_metadata.txt", sep = ',', header=TRUE)
data <- data[!duplicated(genes$gene_name),]
genes <- genes[!duplicated(genes$gene_name),]
row.names(data) <- genes$gene_name


# CASE-CONTROL FEMALE-MALE
c<-samples[samples$condition_name == "naive", ]
cc<-samples[samples$condition_name == "IFNg", ]
control = intersect(c$sample_id, colnames(data))
case= intersect(cc$sample_id, colnames(data))
m<-samples[samples$sex == "male", ]
f<-samples[samples$sex == "female", ]

males = m$sample_id
females = f$sample_id

control_males = intersect(control, males)
control_females = intersect(control, females)

case_males =   intersect(case, males)
case_females = intersect(case, females)


data_NGCs = data[intersect(c(cluster1, cluster2,cluster3), rownames(data)), ]
data_clean <- na.omit(data_NGCs) 
data_clean  = t(data_clean)

data_control_males<- data_clean[rownames(data_clean)    %in% control_males,]
data_control_females<- data_clean[rownames(data_clean)  %in% control_females,]


data_case_males<- data_clean[rownames(data_clean)    %in% case_males,]
data_case_females<- data_clean[rownames(data_clean)  %in% case_females,]


data_cases = rbind(data_case_males,data_case_females )

data_control_males <- data_control_males[,apply(data_control_males, MARGIN = 2, FUN = function(x) sd(x) != 0)]
data_case_males <- data_case_males[,apply(data_case_males, MARGIN = 2, FUN = function(x) sd(x) != 0)]
data_control_females <- data_control_females[,apply(data_control_females, MARGIN = 2, FUN = function(x) sd(x) != 0)]
data_case_females <- data_case_females[,apply(data_case_females, MARGIN = 2, FUN = function(x) sd(x) != 0)]


# control males
c1c <- intersect(cluster1, colnames(data_control_males))
c2c <- intersect(cluster2, colnames(data_control_males))
c3c <- intersect(cluster3, colnames(data_control_males))

cont_males_gene_annotation <- data.frame(clusters = c(rep('cluster1', length(c1c)),rep('cluster2', length(c2c)),rep('cluster3', length(c3c))))

rownames(cont_males_gene_annotation) <- colnames(data_control_males)

# keep Ernest's ordering
cmh = pheatmap(as.data.frame(bicor(data_control_males)), cluster_rows = F, cluster_cols = F, annotation_row = cont_males_gene_annotation, annotation_col = cont_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(cmh,filename="coexp_plots/gaffney_ctrl_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
#cluster based on gaffney
cmhc = pheatmap(as.data.frame(bicor(data_control_males)), cluster_rows = T, cluster_cols = T, annotation_row = cont_males_gene_annotation, annotation_col = cont_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(cmhc,filename="coexp_plots/gaffney_ctrl_males_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#keep Ernest's clusters but change the ordering within the clusters based on gaffney's clustering
cmhc_og <- cmhc$tree_row$labels[cmhc$tree_row$order]
cmhc_og <- c(intersect(cmhc_og, cluster1), intersect(cmhc_og, cluster2), intersect(cmhc_og, cluster3))
cm_gene_annotation <- cont_males_gene_annotation
rownames(cm_gene_annotation) <- cmhc_og
cmhcogdf = as.data.frame(bicor(data_control_males))[cmhc_og,cmhc_og]
cmh_cog = pheatmap(cmhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = cm_gene_annotation, annotation_col = cm_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(cmh_cog,filename="coexp_plots/gaffney_ctrl_males_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)




#control females
c1c <- intersect(cluster1, colnames(data_control_females))
c2c <- intersect(cluster2, colnames(data_control_females))
c3c <- intersect(cluster3, colnames(data_control_females))

cont_females_gene_annotation <- data.frame(clusters = c(rep('cluster1', length(c1c)),rep('cluster2', length(c2c)),rep('cluster3', length(c3c))))

rownames(cont_females_gene_annotation) <- colnames(data_control_females)

# keep Ernest's ordering
cfh = pheatmap(as.data.frame(bicor(data_control_females)), cluster_rows = F, cluster_cols = F, annotation_row = cont_males_gene_annotation, annotation_col = cont_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(cfh,filename="coexp_plots/gaffney_ctrl_females_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
#cluster based on gaffney
cfhc = pheatmap(as.data.frame(bicor(data_control_females)), cluster_rows = T, cluster_cols = T, annotation_row = cont_males_gene_annotation, annotation_col = cont_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(cfhc,filename="coexp_plots/gaffney_ctrl_females_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#keep Ernest's clusters but change the ordering within the clusters based on gaffney's clustering
cfhc_og <- cfhc$tree_row$labels[cfhc$tree_row$order]
cfhc_og <- c(intersect(cfhc_og, cluster1), intersect(cfhc_og, cluster2), intersect(cfhc_og, cluster3))
cf_gene_annotation <- cont_females_gene_annotation
rownames(cf_gene_annotation) <- cfhc_og
cfhcogdf = as.data.frame(bicor(data_control_females))[cfhc_og,cfhc_og]
cfh_cog = pheatmap(cfhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = cf_gene_annotation, annotation_col = cf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(cfh_cog,filename="coexp_plots/gaffney_ctrl_females_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#plot wrt control males

cont_females_wrt_cont_males = as.data.frame(bicor(data_control_females))[cmhc_og,cmhc_og]
cont_females_wrt_cont_males_hm = pheatmap(cont_females_wrt_cont_males, cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(cont_females_wrt_cont_males_hm,filename="coexp_plots/gaffney_control_females_cl_or_based_cont_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#case males

c1c <- intersect(cluster1, colnames(data_case_males))
c2c <- intersect(cluster2, colnames(data_case_males))
c3c <- intersect(cluster3, colnames(data_case_males))

case_males_gene_annotation <- data.frame(clusters = c(rep('cluster1', length(c1c)),rep('cluster2', length(c2c)),rep('cluster3', length(c3c))))

rownames(case_males_gene_annotation) <- colnames(data_case_males)

# keep Ernest's ordering
smh = pheatmap(as.data.frame(bicor(data_case_males)), cluster_rows = F, cluster_cols = F, annotation_row = case_males_gene_annotation, annotation_col = case_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh,filename="coexp_plots/gaffney_case_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
#cluster based on gaffney
smhc = pheatmap(as.data.frame(bicor(data_case_males)), cluster_rows = T, cluster_cols = T, annotation_row = case_males_gene_annotation, annotation_col = case_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(smhc,filename="coexp_plots/gaffney_case_males_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#keep Ernest's clusters but change the ordering within the clusters based on gaffney's clustering
smhc_og <- smhc$tree_row$labels[smhc$tree_row$order]
smhc_og <- c(intersect(smhc_og, cluster1), intersect(smhc_og, cluster2), intersect(smhc_og, cluster3))
cf_gene_annotation <- case_males_gene_annotation
rownames(cf_gene_annotation) <- smhc_og
smhcogdf = as.data.frame(bicor(data_case_males))[smhc_og,smhc_og]
smh_cog = pheatmap(smhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = cf_gene_annotation, annotation_col = cf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh_cog,filename="coexp_plots/gaffney_case_males_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#replot case males based on control males ordering

cf_gene_annotation <- case_males_gene_annotation
cmhc_og <- cmhc$tree_row$labels[cmhc$tree_row$order]
rownames(cf_gene_annotation) <- cmhc_og
smhcogdf = as.data.frame(bicor(data_case_males))[cmhc_og,cmhc_og]
smh_cog = pheatmap(smhcogdf, cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(smh_cog,filename="coexp_plots/gaffney_case_males_cl_or_based_cont_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#case females
c1c <- intersect(cluster1, colnames(data_case_females))
c2c <- intersect(cluster2, colnames(data_case_females))
c3c <- intersect(cluster3, colnames(data_case_females))

case_females_gene_annotation <- data.frame(clusters = c(rep('cluster1', length(c1c)),rep('cluster2', length(c2c)),rep('cluster3', length(c3c))))

rownames(case_females_gene_annotation) <- colnames(data_case_females)

# keep Ernest's ordering
smh = pheatmap(as.data.frame(bicor(data_case_females)), cluster_rows = F, cluster_cols = F, annotation_row = case_females_gene_annotation, annotation_col = case_females_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh,filename="coexp_plots/gaffney_case_females_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
#cluster based on gaffney
smhc = pheatmap(as.data.frame(bicor(data_case_females)), cluster_rows = T, cluster_cols = T, annotation_row = case_females_gene_annotation, annotation_col = case_females_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(smhc,filename="coexp_plots/gaffney_case_females_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


#keep Ernest's clusters but change the ordering within the clusters based on gaffney's clustering
smhc_og <- smhc$tree_row$labels[smhc$tree_row$order]
smhc_og <- c(intersect(smhc_og, cluster1), intersect(smhc_og, cluster2), intersect(smhc_og, cluster3))
cf_gene_annotation <- case_females_gene_annotation
rownames(cf_gene_annotation) <- smhc_og
smhcogdf = as.data.frame(bicor(data_case_females))[smhc_og,smhc_og]
smh_cog = pheatmap(smhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = cf_gene_annotation, annotation_col = cf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh_cog,filename="coexp_plots/gaffney_case_females_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


smhc = pheatmap(as.data.frame(bicor(data_case_males)), cluster_rows = T, cluster_cols = T, annotation_row = case_males_gene_annotation, annotation_col = case_males_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
smhc_og <- smhc$tree_row$labels[smhc$tree_row$order]
smhcogdf = as.data.frame(bicor(data_case_females))[smhc_og,smhc_og]
smh_cog = pheatmap(smhcogdf, cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(smh_cog,filename="coexp_plots/gaffney_case_females_based_case_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)



# replot case females wrt control females
cfhc_og <- cfhc$tree_row$labels[cfhc$tree_row$order]
case_females_wrt_cont_females = as.data.frame(bicor(data_case_females))[cfhc_og,cfhc_og]
case_females_wrt_cont_females_hm = pheatmap(case_females_wrt_cont_females, cluster_rows = F, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(case_females_wrt_cont_females_hm,filename="coexp_plots/gaffney_case_females_cl_or_based_cont_females_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)






smh = pheatmap(as.data.frame(bicor(data_cases_males)), cluster_rows = F, cluster_cols = F, annotation_row = gene_annotation, annotation_col = gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh,filename="coexp_plots/gaffney_case_males_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
smhc = pheatmap(as.data.frame(bicor(data_cases_males)), cluster_rows = T, cluster_cols = T, annotation_row = gene_annotation, annotation_col = gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(smhc,filename="coexp_plots/gaffney_case_males_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)

smhc_og <- smhc$tree_row$labels[smhc$tree_row$order]
smhc_og <- c(intersect(smhc_og, cluster1), intersect(smhc_og, cluster2), intersect(smhc_og, cluster3))
sm_gene_annotation <- gene_annotation
rownames(sm_gene_annotation) <- smhc_og
smhcogdf = as.data.frame(bicor(data_control_males))[smhc_og,smhc_og]
smh_cog = pheatmap(smhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = sm_gene_annotation, annotation_col = sm_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(smh_cog,filename="coexp_plots/gaffney_case_males_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


sfh = pheatmap(as.data.frame(bicor(data_cases_females)), cluster_rows = F, cluster_cols = F, annotation_row = gene_annotation, annotation_col = gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(sfh,filename="coexp_plots/gaffney_case_females_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)
sfhc = pheatmap(as.data.frame(bicor(data_cases_females)), cluster_rows = T, cluster_cols = T, annotation_row = gene_annotation, annotation_col = gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4)
ggsave(sfhc,filename="coexp_plots/gaffney_case_females_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)

sfhc_og <- sfhc$tree_row$labels[sfhc$tree_row$order]
sfhc_og <- c(intersect(sfhc_og, cluster1), intersect(sfhc_og, cluster2), intersect(sfhc_og, cluster3))
sf_gene_annotation <- gene_annotation
rownames(sf_gene_annotation) <- sfhc_og
sfhcogdf = as.data.frame(bicor(data_control_males))[sfhc_og,sfhc_og]
sfh_cog = pheatmap(sfhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = sf_gene_annotation, annotation_col = sf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(sfh_cog,filename="coexp_plots/gaffney_case_females_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


sfmh = pheatmap(as.data.frame(bicor(data_cases)), cluster_rows = T, cluster_cols = T, annotation_row = gene_annotation, annotation_col = gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(sfmh,filename="coexp_plots/gaffney_cases_cl_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)


sfmh_og <- sfmh$tree_row$labels[sfmh$tree_row$order]
sfmh_og <- c(intersect(sfmh_og, cluster1), intersect(sfmh_og, cluster2), intersect(sfmh_og, cluster3))
sf_gene_annotation <- gene_annotation
rownames(sf_gene_annotation) <- sfmh_og
sfmhcogdf = as.data.frame(bicor(data_cases))[sfmh_og,sfmh_og]
sfmh_cog = pheatmap(sfmhcogdf, cluster_rows = F, cluster_cols = F, annotation_row = sf_gene_annotation, annotation_col = sf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(sfmh_cog,filename="coexp_plots/gaffney_cases_cl_or_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)




case_males_cogdf = as.data.frame(bicor(data_cases_males))[sfmh_og,sfmh_og]

case_males = pheatmap(case_males_cogdf, cluster_rows = F, cluster_cols = F, annotation_row = sf_gene_annotation, annotation_col = sf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(males,filename="coexp_plots/gaffney_case_males_cl_order_from_all_cases_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)



case_females_cogdf = as.data.frame(bicor(data_cases_females))[sfmh_og,sfmh_og]

case_females = pheatmap(case_females_cogdf, cluster_rows = F, cluster_cols = F, annotation_row = sf_gene_annotation, annotation_col = sf_gene_annotation, color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_females,filename="coexp_plots/gaffney_case_females_cl_order_from_all_cases_quantile.pdf",width=25,height=25,units="cm",limitsize = FALSE)



diff_cor_cases_females = bicor(data_case_females) - bicor(data_case_males)
diff_cor_cases_males = bicor(data_case_males) - bicor(data_case_females)

diff_cor_controls_females = bicor(data_control_females) - bicor(data_control_males)
diff_cor_controls_males = bicor(data_control_males) - bicor(data_control_females)


case_diff_female_clustered = pheatmap(as.data.frame(diff_cor_cases_females), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_female_clustered,filename="coexp_plots/gaffney_case_diff_females_males_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)


case_diff_male_clustered = pheatmap(as.data.frame(diff_cor_cases_males), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_male_clustered,filename="coexp_plots/gaffney_case_diff_males_females_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)


case_diff_female_clustered = pheatmap(as.data.frame(diff_cor_controls_females), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_female_clustered,filename="coexp_plots/gaffney_control_diff_females_males_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)


case_diff_male_clustered = pheatmap(as.data.frame(diff_cor_controls_males), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_male_clustered,filename="coexp_plots/gaffney_control_diff_males_females_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)


diff_cor_females = bicor(data_case_females) - bicor(data_control_females)
diff_cor_males = bicor(data_case_females) - bicor(data_control_males)


case_diff_female_clustered = pheatmap(as.data.frame(diff_cor_females), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_female_clustered,filename="coexp_plots/gaffney_females_case_cont_diff_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)


case_diff_male_clustered = pheatmap(as.data.frame(diff_cor_males), cluster_rows = T, cluster_cols = T,  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10), fontSize = 4, fontsize_row = 4, fontsize_col = 4, cellheight=4, cellwidth = 4, gaps_row = c(length(c1c), (length(c1c)+length(c2c))), gaps_col = c(length(c1c), (length(c1c)+length(c2c))))
ggsave(case_diff_male_clustered,filename="coexp_plots/gaffney_males_case_cont_diff_cl_quantile_143.pdf",width=25,height=25,units="cm",limitsize = FALSE)



