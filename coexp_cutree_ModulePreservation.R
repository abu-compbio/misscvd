library(pheatmap)
library(ggplot2)
library(WGCNA)
library(DESeq2)

library("RColorBrewer")
setwd('/Users/onurdogan/Desktop/Projects/comp_bio/inflammatory-genes/data/gaffney')

cluster1 <- c('BOC',  'EPHB2',  'NPTXR',  'NDNF',  'NFASC',  'NTN5',  'SEMA3B',  'DSCAM',  'SEMA3A',  'ROBO1',  'DCC',  'EPHA5',  'RELN',  'SEMA3D',  'GFRA1',  'IGSF9',  'PRTG',  'ROBO2',  'SEMA6D',  'RGMA',  'NTRK3',  'FLRT2',  'FZD3',  'BMP7',  'NELL1',  'EPHA10',  'SEMA3E',  'PTK7',  'SEMA4C',  'SEMA4F',  'DRAXIN',  'EPHA1',  'EPHA3',  'CNTN4',  'NPFFR2',  'SEMA6A',  'BDNF',  'EPHA4',  'NTRK2',  'NRG1',  'WNT5A',  'EPHA6',  'CNTN6',  'EPHA7',  'NTNG2',  'PLXNA4',  'UNC5C',  'FGF20',  'FLRT3',  'NRCAM',  'SLIT1',  'EPHB1',  'NTRK1',  'CDH4',  'CDNF',  'SEMA5A',  'UNC5D',  'PDGFRA',  'WNT3',  'GDNF',  'PLXNA3',  'NRP1',  'SHH',  'EPHB4',  'EFNA5',  'LRTM1',  'NTNG1')
cluster2 <- c('GDF15',  'SEMA3C',  'SEMA4B',  'NOTCH2',  'SEMA4D',  'EPHA2',  'UNC5B',  'EDN1',  'EPHB6',  'SEMA3F',  'PLXNB1',  'NOTCH3',  'PDGFRB',  'ALCAM',  'NPTN',  'NENF',  'CXCL12',  'PLXNB2',  'ARRB2',  'PLXNA1',  'CSF1R',  'PLXND1',  'NOTCH1',  'ROBO3',  'BSG',  'LRP1',  'MANF')
cluster3 <- c('FRS2',  'EFNB2',  'SLIT2',  'GDF7',  'NTN4',  'PLXNC1',  'SEMA4A',  'NECTIN1',  'PTPRO',  'ROBO4',  'SEMA5B',  'NGFR',  'PLXNA2',  'NTF3',  'NTN1',  'EFNB1',  'NOG',  'PLXNB3',  'EPHB3',  'NEO1',  'ADAMTSL1',  'SEMA3G',  'LGR6',  'SLIT3',  'ARRB1',  'NRP2',  'GFRA2',  'SCN1B',  'SEMA6B')



data <- read.table('RNA_cqn_matrix.txt')
genes <- read.table('RNA_gene_metadata.txt', header = 1)
samples <- read.table("RNA_sample_metadata.txt", sep = ',', header=TRUE)
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


#analysis starts here
cor <- WGCNA::bicor

reference_data=data_control_males
target_data=data_case_males

corr_control_male=bicor(reference_data)
phetmap_corr_control_male = pheatmap(corr_control_male,cluster_rows = T , cluster_cols = T, fontSize  =2,fontsize_row = 2,fontsize_col =8)
hclust_control_male=phetmap_corr_control_male$tree_row


hclusth = cutree(hclust_control_male, h=4)

module.colours=labels2colors(hclusth)

plotDendroAndColors(hclust_control_male, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='') 

module_eigengenes_control_male = moduleEigengenes(reference_data, module.colours)$eigengenes

multiExpr = list(X = list(data = reference_data), Y = list(data = target_data));
colorList = list(X=module.colours);


mp = modulePreservation(multiExpr,colorList,
                        referenceNetworks = 1,
                        nPermutations = 100,
                        loadPermutedStatistics = FALSE,
                        verbose = 3)
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )  

res=as.data.frame(module.colours)
res['gene']=row.names(corr_control_male)
res['zsummary']=''
Zsummary=as.data.frame(statsZ[, c("Zsummary.pres")])
Zsummary['module']=rownames(statsZ)

for (i in unique(res$module.colours)) {
  res[res$module.colours==i,'zsummary']=Zsummary[Zsummary$module==i,'statsZ[, c("Zsummary.pres")]']
  
}

a=dim(reference_data)[2]

plotMods = !(rownames(statsZ) %in% c("gold"));
text = rownames(statsZ)[plotMods];

colnames(module_eigengenes_control_male)=text
res['eigenCor']=''

#correlation between gene and the corresponding module eigengene

for (i in 1:a){
  gene=res[i,2]
  color=res[i,1]
  eigencor_res=bicor(reference_data[,gene],module_eigengenes_control_male[,color])
  res[i,'eigenCor']=eigencor_res
}

res <- res[order(module.colours),] 

write_xlsx(res, "cutreeCaseFemalevsCaseMale.xlsx")


