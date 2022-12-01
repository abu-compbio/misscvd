library("readxl")
library("writexl")

#Male Female
control_male_vs_case_male <- read_excel("gaffney_dist_filtered_case_males_ctrl_males.xlsx")
control_female_vs_case_female <- read_excel("gaffney_dist_filtered_case_females_ctrl_females.xlsx")

control_male_vs_case_male=as.data.frame(control_male_vs_case_male)
control_female_vs_case_female=as.data.frame(control_female_vs_case_female)

control_male_vs_case_male$RankMale<-rank(control_male_vs_case_male$NA_col)
control_female_vs_case_female$RankFemale<-rank(control_female_vs_case_female$NA_col)

s <- merge(control_male_vs_case_male,control_female_vs_case_female,by = 'gene', all = TRUE)

s$diff=abs(s$RankMale-s$RankFemale)
res=s[order(s$diff, decreasing = TRUE),]

write_xlsx(data.frame(res), "male_female-most-diff_genes.xlsx")

# Control Case
control_male_vs_control_female <- read_excel("gaffney_dist_filtered_ctrl_males_ctrl_females.xlsx")
case_male_vs_case_female <- read_excel("gaffney_dist_filtered_case_males_case_females.xlsx")

control_male_vs_control_female=as.data.frame(control_male_vs_control_female)
case_male_vs_case_female= as.data.frame(case_male_vs_case_female)

control_male_vs_control_female$RankControl<-rank(control_male_vs_control_female$NA_col)
case_male_vs_case_female$RankFCase<-rank(case_male_vs_case_female$NA_col)

s <- merge(control_male_vs_control_female,case_male_vs_case_female,by = 'gene', all = TRUE)

s$diff=abs(s$RankControl-s$RankFCase)

res=s[order(s$diff, decreasing = TRUE),]

write_xlsx(data.frame(res), "control_case-most-diff_genes.xlsx")

