install.packages("tidyverse")
install.packages("dplyr")
devtools::install_github("MRCIEU/TwoSampleMR")

library(data.table)
library(tidyverse)
library(dplyr)
library(TwoSampleMR)
library(reshape2)


##########################################################################################
##FUNCTIONS
##########################################################################################

is_outlier <- function(x, iqrfac = 3) {
  quants <- quantile(x, na.rm = TRUE)
  iqr <- quants[4] - quants[2]
  !is.na(x) & (x < (quants[2] - iqrfac*iqr) | (quants[4] + iqrfac*iqr) < x)
}

##########################################################################################
#Set directory where phenotypic data is stored
##########################################################################################

setwd("~/PRS_AD_BRAIN_STRUCTURES/prs_alzheimers_brain_structures/Data/UKBIOBANK/Phenotypic_data/")
cortical_structures_area <- fread("cortical_area_ukbiobank192.txt")

##########################################################################################
## return indices of repeat visits
##########################################################################################

ind <- cortical_structures_area[,grep(".3.0",colnames(cortical_structures_area))]
160 162 164 166 168 170 172 174 176 178 180

##########################################################################################
##Remove repeat image data
##########################################################################################

cortical_structures_area <- cortical_structures_area[,-c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,139,140,142,154:185,187)]
colnames(cortical_structures_area)[3:83] <- c("l_total_area","r_total_area","l_banksst_area","r_banksst_area","l_caudalanteriorcingulate_area","r_caudalanteriorcingulate_area","l_caudalmiddlefrontal_area","r_caudalmiddlefrontal_area","l_cuneus_area","r_cuneus_area","l_entorhinal_area","r_entorhinal_area","l_frontalpole_area","r_frontalpole_area","l_fusiform_area","r_fusiform_area","l_inferiorparietal_area","r_inferiorparietal_area","l_inferiortemporal_area","r_inferiortemporal_area","l_insula_area","r_insula_area","l_isthmuscingulate_area","r_isthmuscingulate_area","l_lateraloccipital_area","r_lateraloccipital_area","l_lateralorbitofrontal_area","r_lateralorbitofrontal_area","l_lingual_area","r_lingual_area","l_medialorbitofrontal_area","r_medialorbitofrontal_area","l_middletemporal_area","r_middletemporal_area","l_paracentral_area","r_paracentral_area","l_parahippocampal_area","r_parahippocampal_area","l_parsopercularis_area","r_parsopercularis_area","l_parsorbitalis_area","r_parsorbitalis_area","l_parstriangularis_area","r_parstriangularis_area","l_pericalcarine_area","r_pericalcarine_area","l_postcentral_area","r_postcentral_area","l_posteriorcingulate_area","r_posteriorcingulate_area","l_precentral_area","r_precentral_area","l_precuneus_area","r_precuneus_area","l_rostralanteriorcingulate_area","r_rostralanteriorcingulate_area","l_rostralmiddlefrontal_area","r_rostralmiddlefrontal_area","l_superiorfrontal_area","r_superiorfrontal_area","l_superiorparietal_area","r_superiorparietal_area","l_superiortemporal_area","r_superiortemporal_area","l_supramarginal_area","r_supramarginal_area","l_transversetemporal_area","r_transversetemporal_area","age_imaging","sex","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","assessment_centre")

##########################################################################################
##Merge in genetic data
##########################################################################################

intelligence_SNPS_UKB <- fread("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/ukbb_intelligence_snps_proxies_hill.raw")
colnames(intelligence_SNPS_UKB)[2] <- "appieu"
linker_file <- fread("~/48970/dev/release_candidate/data/linker_file.csv")
intelligence_SNPS_UKB_phenid <-  merge(linker_file,intelligence_SNPS_UKB,by="appieu")
colnames(intelligence_SNPS_UKB_phenid)[2] <- "eid"
cortical_area_gen <- merge(intelligence_SNPS_UKB_phenid,cortical_structures_area,by.x="eid",by.y="projectID")
names(cortical_area_gen) <- gsub("\\(..)", "",names(cortical_area_gen))

##########################################################################################
##Remove outliers if they lie below or above 3*IQR
##########################################################################################

cortical_area_gen[,c(162:229)] <- Map(replace, cortical_area_gen[,c(162:229)], lapply(cortical_area_gen[,c(162:229)], is_outlier), NA)

cortical_area_complete<- cortical_area_gen[apply(!is.na(cortical_area_gen[,162:229]),1,any),]

##########################################################################################
#Order data by age and stratify into three equal sub-samples
##########################################################################################

cortical_area_complete <- cortical_area_complete[order(cortical_area_complete$age_imaging),]
cortical_area_completeT1 <-cortical_area_complete[1:9377,]
cortical_area_completeT2 <- cortical_area_complete[9378:18753,]
cortical_area_completeT3 <- cortical_area_complete[18754:28129,]

##########################################################################################
#Select brain imaging columns
##########################################################################################

cortical_area_subset1 <- cortical_area_completeT1[,c(162:229)]
cortical_area_subset2 <- cortical_area_completeT2[,c(162:229)]
cortical_area_subset3 <- cortical_area_completeT3[,c(162:229)]


##########################################################################################
##calculate mean for cortical area structures [age tertile 1]
##########################################################################################

out_cort_area1<-t(cortical_area_subset1) %>%
data.frame() %>%
group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
summarise_all(mean) %>%
data.frame() %>%
column_to_rownames(var ='id') %>%
t()

##########################################################################################
##calculate mean for cortical area structures [age tertile 2]
##########################################################################################

out_cort_area2<- t(cortical_area_subset2) %>%
data.frame() %>%
group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
summarise_all(mean) %>%
data.frame() %>%
column_to_rownames(var ='id') %>%
t()

##########################################################################################
##calculate mean for cortical area structures [age tertile 3]
##########################################################################################

out_cort_area3 <-t(cortical_area_subset3) %>%
data.frame() %>%
group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
summarise_all(mean) %>%
data.frame() %>%
column_to_rownames(var ='id') %>%
t()

##########################################################################################
##standardise mean measures
##########################################################################################

data.standardised_cortical_area1 <-as.data.frame(scale(out_cort_area1))
names(data.standardised_cortical_area1) <- sapply(names(data.standardised_cortical_area1), function(x) {
  y <- paste0(x, "_standardised_T1")
}
)

##########################################################################################
##standardise mean measures
##########################################################################################

data.standardised_cortical_area2 <-as.data.frame(scale(out_cort_area2))
names(data.standardised_cortical_area2) <- sapply(names(data.standardised_cortical_area2), function(x) {
  y <- paste0(x, "_standardised_T2")
}
)

##########################################################################################
##standardise mean measures
##########################################################################################

data.standardised_cortical_area3 <-as.data.frame(scale(out_cort_area3))
names(data.standardised_cortical_area3) <- sapply(names(data.standardised_cortical_area3), function(x) {
  y <- paste0(x, "_standardised_T3")
}
)

##########################################################################################
#Bind genetic and phenotypic data
##########################################################################################

genetic_standardised_cortical_area1 <- cbind(cortical_area_completeT1[,c(1,8:161,230:241)],data.standardised_cortical_area1)
genetic_standardised_cortical_area2 <- cbind(cortical_area_completeT2[,c(1,8:161,230:241)],data.standardised_cortical_area2)
genetic_standardised_cortical_area3 <- cbind(cortical_area_completeT3[,c(1,8:161,230:241)],data.standardised_cortical_area3)

##########################################################################################
##filter dataset to include complete data to be used in regression models - Tertile 1
##########################################################################################

genetic_cortical_area1 <- genetic_standardised_cortical_area1[!is.na(genetic_standardised_cortical_area1$sex) & !is.na(genetic_standardised_cortical_area1$age_imaging) & !is.na(genetic_standardised_cortical_area1$PC1)]
descriptive_stats_cortical_area1 <- pairwiseCount(genetic_cortical_area1[,c(2:154)],genetic_cortical_area1[,c(168:201)])
descriptive_stats_cortical_area1 <- t(descriptive_stats_cortical_area1)
descriptive_stats_cortical_area_reshaped1<-setNames(melt(descriptive_stats_cortical_area1),c('Outcome','SNP','N'))
write.table(descriptive_stats_cortical_area_reshaped1,"~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T1_intelligence_snps.csv",row.names=FALSE, quote=FALSE,sep=",")

##########################################################################################
##filter dataset to include complete data to be used in regression models - Tertile 2
##########################################################################################

genetic_cortical_area2 <- genetic_standardised_cortical_area2[!is.na(genetic_standardised_cortical_area2$sex) & !is.na(genetic_standardised_cortical_area2$age_imaging) & !is.na(genetic_standardised_cortical_area2$PC1)]
descriptive_stats_cortical_area2 <- pairwiseCount(genetic_cortical_area2[,c(2:154)],genetic_cortical_area2[,c(168:201)])
descriptive_stats_cortical_area2 <- t(descriptive_stats_cortical_area2)
descriptive_stats_cortical_area_reshaped2<-setNames(melt(descriptive_stats_cortical_area2),c('Outcome','SNP','N'))
write.table(descriptive_stats_cortical_area_reshaped2,"~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T2_intelligence_snps.csv",row.names=FALSE, quote=FALSE,sep=",")

##########################################################################################
##filter dataset to include complete data to be used in regression models - Tertile 3
##########################################################################################

genetic_cortical_area3 <- genetic_standardised_cortical_area3[!is.na(genetic_standardised_cortical_area3$sex) & !is.na(genetic_standardised_cortical_area3$age_imaging) & !is.na(genetic_standardised_cortical_area3$PC1)]
descriptive_stats_cortical_area3 <- pairwiseCount(genetic_cortical_area3[,c(2:154)],genetic_cortical_area3[,c(168:201)])
descriptive_stats_cortical_area3 <- t(descriptive_stats_cortical_area3)
descriptive_stats_cortical_area_reshaped3 <-setNames(melt(descriptive_stats_cortical_area3),c('Outcome','SNP','N'))
write.table(descriptive_stats_cortical_area_reshaped3,"~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T3_intelligence_snps.csv",row.names=FALSE, quote=FALSE,sep=",")

##########################################################################################
#SET RESULTS DIRECTORY
#############################################################################################################################################

setwd("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/GWAS")

##########################################################################################
##Select genetic variants for ivs_vec
##########################################################################################

ivs_vec <- names(genetic_standardised_cortical_area1)[2:154]

##########################################################################################
##Select brain imaging measures for dvs_vec
##########################################################################################

dvs_vec <- names(genetic_standardised_cortical_area1)[168:201]

##########################################################################################
##Select covariates for covs_vec
##########################################################################################

covs_vec <- paste("age_imaging","+","total_area_standardised_T1","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5","+","PC6","+","PC7","+","PC8","+","PC9","+","PC10")

#############################################################################################################################################
##RUN REGRESSION MODELS
#############################################################################################################################################
ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=genetic_standardised_cortical_area1)
})

##########################################################################################
# Creating / combining results
##########################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)

##########################################################################################
# subset results just for MR
##########################################################################################
x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=FALSE); 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
sample_size <- read.csv("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T1_intelligence_snps.csv")

intelligence_alleles <- read.table("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/Discovery_samples/SNPs/Intelligence/Proxies/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
mr_dataset <- merge(combined_results,sample_size,by=c("SNP","Outcome"))
mr_dataset$SNP <- sub("_[^_]+$", "",mr_dataset$SNP)
mr_dataset <- merge(mr_dataset,intelligence_alleles,by="SNP")
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1],".txt"),quote = FALSE,row.names = FALSE))

########################################################################################################################################################
#TERTILE 2
##########################################################################################

##########################################################################################
##Select genetic variants for ivs_vec
##########################################################################################

ivs_vec <- names(genetic_standardised_cortical_area2)[2:154]

##########################################################################################
##Select brain imaging measures for dvs_vec
##########################################################################################
dvs_vec <- names(genetic_standardised_cortical_area2)[168:201]

##########################################################################################
##Select columns for covariates
##########################################################################################
covs_vec <- paste("age_imaging","+","sex","+","total_area_standardised_T2","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5","+","PC6","+","PC7","+","PC8","+","PC9","+","PC10")
ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)


##########################################################################################
##RUN REGRESSION MODELS
##########################################################################################

lm_results <- lapply(formulas, function(x) {
  lm(x, data=genetic_standardised_cortical_area2)
})


##########################################################################################
# Creating / combining results
##########################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)

###########################################################################################
#subset results just for MR
##########################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=FALSE); 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
sample_size <- read.csv("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T2_intelligence_snps.csv")

intelligence_alleles <- read.table("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/Discovery_samples/SNPs/Intelligence/Proxies/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
mr_dataset <- merge(combined_results,sample_size,by=c("SNP","Outcome"))
mr_dataset$SNP <- sub("_[^_]+$", "",mr_dataset$SNP)
mr_dataset <- merge(mr_dataset,intelligence_alleles,by="SNP")
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1],".txt"),quote = FALSE,row.names = FALSE))

####################################################################################################################################################
#TERTILE 3
####################################################################################################################################################

##########################################################################################
##Select genetic variants for ivs_vec
##########################################################################################

ivs_vec <- names(genetic_standardised_cortical_area2)[2:154]

##########################################################################################
##Select brain imaging measures for dvs_vec
##########################################################################################
dvs_vec <- names(genetic_standardised_cortical_area2)[168:201]

covs_vec <- paste("age_imaging","+","sex","+","total_area_standardised_T3","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5","+","PC6","+","PC7","+","PC8","+","PC9","+","PC10")
ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

##########################################################################################
#RUN REGRESSION MODELS
##########################################################################################

lm_results <- lapply(formulas, function(x) {
  lm(x, data=genetic_standardised_cortical_area3)
})

##########################################################################################
# Creating / combining results
##########################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)

##########################################################################################
# subset results just for MR
##########################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=FALSE); 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
sample_size <- read.csv("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/UK_Biobank/GeneticData/Intelligence/Descriptive_statistics/descriptive_stats_cortical_area_ukbiobank_T3_intelligence_snps.csv")

intelligence_alleles <- read.table("~/MR_EDUCATION_IQ_BRAIN/prs_edu_iq_brain/Data/Discovery_samples/SNPs/Intelligence/Proxies/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
mr_dataset <- merge(combined_results,sample_size,by=c("SNP","Outcome"))
mr_dataset$SNP <- sub("_[^_]+$", "",mr_dataset$SNP)
mr_dataset <- merge(mr_dataset,intelligence_alleles,by="SNP")
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1],".txt"),quote = FALSE,row.names = FALSE))

