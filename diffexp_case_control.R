setwd('/Users/onurdogan/Desktop/diffexp-results_analysis/Gaffney/')
library("writexl")

#case-control
mf_case <- read.csv("DESeq2_ernest_IFNg.csv")
mf_control <- read.csv("DESeq2_ernest_naive.csv")

rownames(mf_case)=mf_case$X
rownames(mf_control)=mf_control$X

mf_case_05<-mf_case[which(mf_case[,'padj']<0.05,),'X']
mf_control_05<-mf_control[which(mf_control[,'padj']<0.05,),'X']

diff_pval.05<-setdiff(mf_case_05,mf_control_05)

res_mf_case.05=mf_case[diff_pval.05,c('X','log2FoldChange','padj')]
res_mf_control.05=mf_control[diff_pval.05,c('log2FoldChange','padj')]

colnames(res_mf_case.05)=c("X","fcase_vs_mcase_LFC","fcase_vs_mcase_padj")
colnames(res_mf_control.05)=c("fcontrol_vs_mcontrol_LFC","fcontrol_vs_mcontrol_padj")

s=cbind(res_mf_case.05,res_mf_control.05)

write_xlsx(data.frame(s), "Gaffney-case_control-diff-0.05.xlsx")

#.1
mf_case_1<-mf_case[which(mf_case[,'padj']<0.1,),'X']
mf_control_1<-mf_control[which(mf_control[,'padj']<0.1,),'X']

diff_pval.1<-setdiff(mf_case_1,mf_control_1)

res_mf_case.1=mf_case[diff_pval.1,c('X','log2FoldChange','padj')]
res_mf_control.1=mf_control[diff_pval.1,c('log2FoldChange','padj')]

colnames(res_mf_case.1)=c("X","fcase_vs_mcase_LFC","fcase_vs_mcase_padj")
colnames(res_mf_control.1)=c("fcontrol_vs_mcontrol_LFC","fcontrol_vs_mcontrol_padj")

s2=cbind(res_mf_case.1,res_mf_control.1)

write_xlsx(data.frame(s2), "Gaffney-case_control-diff-0.1.xlsx")

#control-case

mf_case_05<-mf_case[which(mf_case[,'padj']<0.05,),'X']
mf_control_05<-mf_control[which(mf_control[,'padj']<0.05,),'X']

diff_pval.05<-setdiff(mf_control_05,mf_case_05)

res_mf_case.05=mf_case[diff_pval.05,c('log2FoldChange','padj')]
res_mf_control.05=mf_control[diff_pval.05,c('X','log2FoldChange','padj')]

colnames(res_mf_case.05)=c("fcase_vs_mcase_LFC","fcase_vs_mcase_padj")
colnames(res_mf_control.05)=c("X","fcontrol_vs_mcontrol_LFC","fcontrol_vs_mcontrol_padj")

s=cbind(res_mf_control.05,res_mf_case.05)

write_xlsx(data.frame(s), "Gaffney-control_case-diff-0.05.xlsx")

#.1
mf_case_1<-mf_case[which(mf_case[,'padj']<0.1,),'X']
mf_control_1<-mf_control[which(mf_control[,'padj']<0.1,),'X']

diff_pval.1<-setdiff(mf_control_1,mf_case_1)

res_mf_case.1=mf_case[diff_pval.1,c('log2FoldChange','padj')]
res_mf_control.1=mf_control[diff_pval.1,c('X','log2FoldChange','padj')]

colnames(res_mf_case.1)=c("fcase_vs_mcase_LFC","fcase_vs_mcase_padj")
colnames(res_mf_control.1)=c("X","fcontrol_vs_mcontrol_LFC","fcontrol_vs_mcontrol_padj")

s2=cbind(res_mf_control.1,res_mf_case.1)

write_xlsx(data.frame(s2), "Gaffney-control_case-diff--0.1.xlsx")



