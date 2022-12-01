setwd('/Users/onurdogan/Desktop/diffexp-results_analysis/Gaffney/')
library("writexl")

#diffexp analysis part
# male-female
male <- read.csv("DESeq2_ernest_male.csv")
female <- read.csv("DESeq2_ernest_female.txt")

rownames(male)=male$X
rownames(female)=female$X

cc_male_05<-male[which(male[,'padj']<0.05,),'X']
cc_female_05<-female[which(female[,'padj']<0.05,),'X']

diff_pval.05<-setdiff(cc_male_05,cc_female_05)

res_cc_male.05=male[diff_pval.05,c('X','log2FoldChange','padj')]
res_cc_female.05=female[diff_pval.05,c('log2FoldChange','padj')]

colnames(res_cc_male.05)=c("X","mcase_vs_mcontrol_LFC","mcase_vs_mcontrol_padj")
colnames(res_cc_female.05)=c("fcase_vs_fcontrol_LFC","fcase_vs_fcontrol_padj")

s=cbind(res_cc_male.05,res_cc_female.05)

write_xlsx(data.frame(s), "Gaffney-male_female-diff-0.05.xlsx")

#.1
cc_male_1<-male[which(male[,'padj']<0.1,),'X']
cc_female_1<-female[which(female[,'padj']<0.1,),'X']

diff_pval.1<-setdiff(cc_male_1,cc_female_1)

res_cc_male.1=male[diff_pval.1,c('X','log2FoldChange','padj')]
res_cc_female.1=female[diff_pval.1,c('log2FoldChange','padj')]

colnames(res_cc_male.1)=c("X","mcase_vs_mcontrol_LFC","mcase_vs_mcontrol_padj")
colnames(res_cc_female.1)=c("fcase_vs_fcontrol_LFC","fcase_vs_fcontrol_padj")

s2=cbind(res_cc_male.1,res_cc_female.1)

write_xlsx(data.frame(s2), "Gaffney-male_female-diff-0.1.xlsx")


# sign analysis part


male_05<-male[which(male[,'padj']<0.05,),'X']
female_05<-female[which(female[,'padj']<0.05,),'X']

intersect_05<-intersect(male_05,female_05)

a=male[intersect_05,c('X','log2FoldChange')]
b=female[intersect_05,c('X','log2FoldChange')]

for (i in intersect_05){
  if(sign(a[i,'log2FoldChange'])!=sign(b[i,'log2FoldChange'])){
    print(a[i,'X'])
  }
}

male_1<-male[which(male[,'padj']<0.1,),'X']
female_1<-female[which(female[,'padj']<0.1,),'X']

intersect_1<-intersect(male_1,female_1)

a=male[intersect_1,c('X','log2FoldChange')]
b=female[intersect_1,c('X','log2FoldChange')]

for (i in intersect_1){
  if(sign(a[i,'log2FoldChange'])!=sign(b[i,'log2FoldChange'])){
    print(a[i,'X'])
  }
}



#female-male 
# male-female

cc_male_05<-male[which(male[,'padj']<0.05,),'X']
cc_female_05<-female[which(female[,'padj']<0.05,),'X']

diff_pval.05<-setdiff(cc_female_05,cc_male_05)

res_cc_male.05=male[diff_pval.05,c('log2FoldChange','padj')]
res_cc_female.05=female[diff_pval.05,c('X','log2FoldChange','padj')]

colnames(res_cc_male.05)=c("mcase_vs_mcontrol_LFC","mcase_vs_mcontrol_padj")
colnames(res_cc_female.05)=c("X","fcase_vs_fcontrol_LFC","fcase_vs_fcontrol_padj")

s=cbind(res_cc_female.05,res_cc_male.05)

write_xlsx(data.frame(s), "Gaffney-female_male-diff-0.05.xlsx")

#.1
cc_male_1<-male[which(male[,'padj']<0.1,),'X']
cc_female_1<-female[which(female[,'padj']<0.1,),'X']

diff_pval.1<-setdiff(cc_female_1,cc_male_1)

res_cc_male.1=male[diff_pval.1,c('log2FoldChange','padj')]
res_cc_female.1=female[diff_pval.1,c('X','log2FoldChange','padj')]

colnames(res_cc_male.1)=c("mcase_vs_mcontrol_LFC","mcase_vs_mcontrol_padj")
colnames(res_cc_female.1)=c("X","fcase_vs_fcontrol_LFC","fcase_vs_fcontrol_padj")

s2=cbind(res_cc_female.1,res_cc_male.1)

write_xlsx(data.frame(s2), "Gaffney-female_male-diff-0.1.xlsx")


# sign analysis part

male_05<-male[which(male[,'padj']<0.05,),'X']
female_05<-female[which(female[,'padj']<0.05,),'X']

intersect_05<-intersect(female_05,male_05)

a=male[intersect_05,c('X','log2FoldChange')]
b=female[intersect_05,c('X','log2FoldChange')]

for (i in intersect_05){
  if(sign(a[i,'log2FoldChange'])!=sign(b[i,'log2FoldChange'])){
    print(b[i,'X'])
  }
}

male_1<-male[which(male[,'padj']<0.1,),'X']
female_1<-female[which(female[,'padj']<0.1,),'X']

intersect_1<-intersect(female_1,male_1)

a=male[intersect_1,c('X','log2FoldChange')]
b=female[intersect_1,c('X','log2FoldChange')]

for (i in intersect_1){
  if(sign(a[i,'log2FoldChange'])!=sign(b[i,'log2FoldChange'])){
    print(b[i,'X'])
  }
}
