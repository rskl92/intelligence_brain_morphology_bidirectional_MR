library(meta)
library(readxl)
setwd("~/PRS_EDU_COG/Results/Intelligence/Meta_analysis_adulthood")
meta_analysis_data <- read_excel("mr_intelligence_metanalysis_adulthood_0305022_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Inverse variance weighted",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_ivw_adulthood_intelligence.txt")
print(metanalysis)
sink()

meta_analysis_data <- read_excel("mr_intelligence_metanalysis_adulthood_0305022_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="MR Egger",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_mr_egger_adulthood_intelligence.txt")
print(metanalysis)
sink()

meta_analysis_data <- read_excel("mr_intelligence_metanalysis_adulthood_0305022_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Weighted median",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_mrweightedmedian_adulthood_intelligence.txt")
print(metanalysis)
sink() 

meta_analysis_data <- read_excel("mr_intelligence_metanalysis_adulthood_0305022_updated.xlsx")
meta_analysis_data <- meta_analysis_data[meta_analysis_data$method=="Weighted mode",]
metaanalysis <- metagen(TE,seTE,data=meta_analysis_data,studlab = Study,comb.random=TRUE, method.tau="DL",prediction = TRUE,sm="SMD")
metanalysis <- update.meta(metaanalysis,byvar=meta_analysis_data$outcome,comb.random=TRUE,comb.fixed=FALSE)
sink("results_metaanalysis_weightedmode_adulthood_intelligence.txt")
print(metanalysis)
sink()