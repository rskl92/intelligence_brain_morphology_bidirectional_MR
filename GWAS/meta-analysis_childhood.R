library(meta)
library(readxl)
setwd("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_childhood")
meta_analysis_data <- read_excel("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_childhood/mr_results_intelligence_brain_030522_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Inverse variance weighted",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_ivw_childhood_intelligence.txt")
print(metanalysis)
sink()
pval <- as.data.frame(metanalysis$pval.random.w)
write.csv(pval,"pvalues_metanalysis.csv",row.names=FALSE,quote=FALSE)

meta_analysis_data <- read_excel("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_childhood/mr_results_intelligence_brain_030522_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="MR Egger",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_mr_egger_childhood_intelligence.txt")
print(metanalysis)
sink()

meta_analysis_data <- read_excel("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_childhood/mr_results_intelligence_brain_030522_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Weighted median",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_mrweightedmedian_childhood_intelligence.txt")
print(metanalysis)
sink() 

meta_analysis_data <- read_excel("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_childhood/mr_results_intelligence_brain_030522_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Weighted mode",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$Outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_weightedmode_childhood_intelligence.txt")
print(metanalysis)
sink()