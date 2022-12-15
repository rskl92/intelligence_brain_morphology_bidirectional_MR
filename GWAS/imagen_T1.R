#Load packages
library(haven)
library(repurrrsive)
library(purrr)
library(stringr)
library(broom)/
library(broom)
library(data.table)
library(readr)
library(tidyverse)
library(readxl)
library(widyr)
library(psych)

##Importing IMAGEN genetic data
IMAGEN_intelligence_snps <- fread("~/PRS_EDU_COG/Data/IMAGEN/Genetic/Intelligence/SNPs/imagen_intelligence_hill.txt")
IMAGEN_frequencies <- fread("~/PRS_EDU_COG/Data/IMAGEN/Genetic/Intelligence/SNPs/IMAGEN_intelligence_frequencies.txt")
principal_components <- fread("~/IMAGENQC/imagen.merge.all_phase3.no_duplicate_variants.clean.no_ac_gt_snps.eigenvec")
intelligence_alleles <- read.table("~/PRS_EDU_COG/Data/Discovery_Samples/Intelligence/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=T)
five_principal_components <- principal_components[,c(1,3:7)]
five_principal_components$V1 <- gsub("IMAGEN","",five_principal_components$V1)
five_principal_components$V1 <- as.numeric(five_principal_components$V1)
colnames(five_principal_components) <- c("ID","PC1","PC2","PC3","PC4","PC5")
colnames(IMAGEN_intelligence_snps)[1] <- "ID"
five_principal_components <- five_principal_components[2505:4339,]
failed_ancestry_check <- fread("~/IMAGENQC/non_europeans_imagen.txt")
genetic_data<- merge(IMAGEN_intelligence_snps,five_principal_components,by="ID")
colnames(genetic_data)[1]<- "ID"
IMAGEN_siblings <- read.table("~/IMAGEN/freesurfer_processed/BL/IMAGEN_siblings.txt",he=F)
genetic_data <- genetic_data[!genetic_data$ID %in% failed_ancestry_check$FID,]
genetic_data <- genetic_data[!genetic_data$ID %in% IMAGEN_siblings$V1,]

##SET WORKING DIRECTORY
setwd("~/IMAGEN/freesurfer_processed/BL")
list.files()

##IMPORTING IMAGEN MEASURES FROM MRI DATA
thickness_left_hemisphere <- read.table("lh.aparc.thickness.tsv",he=T)
colnames(thickness_left_hemisphere)[1] <- "ID"
thickness_right_hemisphere <- read.table("rh.aparc.thickness.tsv",he=T)
colnames(thickness_right_hemisphere)[1] <- "ID"
area_left_hemisphere <- read.delim("lh.aparc.area.tsv",he=T)
colnames(area_left_hemisphere)[1] <- "ID"
area_right_hemisphere <- read.table("rh.aparc.area.tsv",he=T)
colnames(area_right_hemisphere)[1] <- "ID"

##Demographic variables
age_months<- fread("~/IMAGEN/freesurfer_processed/BL/IMAGEN_age_at_mri_acquisition_BL.csv")
age_months <- age_months[,c(1,5)]
colnames(age_months)[1:2] <- c("ID","age")
age_months$age <- as.numeric(age_months$age)
sex <- read.csv("~/IMAGEN/freesurfer_processed/BL/id_sex.csv")
sex$gender2[sex$gender=="male"] <-0
sex$gender2[sex$gender=="female"] <-1
sex$gender2[sex$gender=="unknown"] <- NA
colnames(sex)[3] <- "sex"
sex <- sex[,c(1,3)]
sex$sex <- as.factor(sex$sex)

##Freesurfer vars
freesurfer_status <- read.table("~/IMAGEN/freesurfer_processed/BL/freesurfer_status_edited.txt",he=T) ##image acquisition failures
euler_metric <- read.table("euler.tsv",he=T)
subcortical_structures <- read.table("~/IMAGEN/freesurfer_processed/BL/aseg.volume.tsv",he=T)
colnames(subcortical_structures)[1] <- "ID" 
##merging MRI datasets
left_hemisphere <- merge(area_left_hemisphere,thickness_left_hemisphere,by="ID")
right_hemisphere <- merge(area_right_hemisphere,thickness_right_hemisphere,by="ID")
cortical_structures <- merge(left_hemisphere,right_hemisphere,by="ID")

##Remove IMAGEN phenotypes with value of 0
is.na(cortical_structures) <- !cortical_structures
is.na(subcortical_structures) <- !subcortical_structures

##Merge with failure acquisition status
cortical_structures <- merge(cortical_structures,freesurfer_status,by="ID")
subcortical_structures <- merge(subcortical_structures,freesurfer_status,by="ID")


##Keep images succeeding acquisition
cortical_structures <- cortical_structures[cortical_structures$Exitcode==0,]
subcortical_structures<- subcortical_structures[subcortical_structures$Exitcode==0,]
############################################################################################################################################

##Remove outliers if they lie below or above 3*IQR
##WRITE FUNCTION
is_outlier <- function(x, iqrfac = 3) {
  quants <- quantile(x, na.rm = TRUE)
  iqr <- quants[4] - quants[2]
  !is.na(x) & (x < (quants[2] - iqrfac*iqr) | (quants[4] + iqrfac*iqr) < x)
}

lapply(cortical_structures[,c(2:141)],is_outlier)
cortical_structures[,c(2:141)] <- Map(replace, cortical_structures[,c(2:141)], lapply(cortical_structures[,c(2:141)], is_outlier), NA)
lapply(subcortical_structures[,c(2,6:9,12:14,16,20,24:30,53:55,67)],is_outlier)
subcortical_structures[,c(2,6:9,12:14,16,20,23:30,53:55,67)]<- Map(replace, subcortical_structures[,c(2,6:9,12:14,16,20,23:30,53:55,67)], lapply(subcortical_structures[,c(2,6:9,12:14,16,20,23:30,53:55,67)], is_outlier), NA)



##Merge genetic and imaging datasets 
cortical_structures_snps <- merge(cortical_structures,age_months,by="ID")
cortical_structures_snps <- merge(cortical_structures_snps,sex,by="ID")
cortical_structures_snps <- merge(genetic_data,cortical_structures_snps,by="ID")


##Merge genetic and imaging datasets 
subcortical_structures_snps <- merge(genetic_data,subcortical_structures,by="ID")
subcortical_structures_snps <- merge(subcortical_structures_snps,age_months,by="ID")
subcortical_structures_snps <- merge(subcortical_structures_snps,sex,by="ID")



##select ROIs
cortical_structures2<- cortical_structures_snps[,c(160:193,195:263,265:299)]
total_surface_area <- rowSums(cortical_structures_snps[,c(160:193,230:263)])
subcortical_structures2<- subcortical_structures_snps[,c(160,164:167,170:172,174,178,182:188,211:212,225)]

##rename white matter
colnames(subcortical_structures2)[18] <- "l_CorticalWhiteMatterVol"  
colnames(subcortical_structures2)[19] <- "r_CorticalWhiteMatterVol"  




##calculate mean for cortical structures
out_cort<-t(cortical_structures2) %>%
  data.frame() %>%
  group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
  summarise_all(mean) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t()



##calculate mean for subcortical structures
out_subcort<-t(subcortical_structures2) %>%
  data.frame() %>%
  group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
  summarise_all(mean) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t()

out_cort2 <- cbind(out_cort,total_surface_area)

##########################################################################################################################################
##standardise mean measures
##########################################################################################################################################

##standardise mean measures - cortical
data.standardised_cortical <- scale(out_cort2)
data.standardised_cortical<- as.data.frame(data.standardised_cortical)
setnames(data.standardised_cortical,paste0(names(data.standardised_cortical),"_standardised"))

##standardise mean measures - subcortical
data.standardised_subcortical <- scale(out_subcort)
data.standardised_subcortical<- as.data.frame(data.standardised_subcortical)
setnames(data.standardised_subcortical,paste0(names(data.standardised_subcortical),"_standardised"))

##########################################################################################################################################
#bind standardised measures with SNPs and covariates
##########################################################################################################################################

cortical_structures_snps2<- cbind(cortical_structures_snps[,c(1:159,301:302)], data.standardised_cortical)
subcortical_structures_snps2<- cbind(subcortical_structures_snps[,c(1:159,227:228)], data.standardised_subcortical)

#########################################################################################################################################
##remove (/ from rsids)
#########################################################################################################################################
names(cortical_structures_snps2) <- gsub("\\(..)", "", names(cortical_structures_snps2))
names(subcortical_structures_snps2) <- gsub("\\(..)", "", names(subcortical_structures_snps2))
colnames(cortical_structures_snps2)[160] <- "age_imaging"
colnames(subcortical_structures_snps2)[160] <- "age_imaging"
#########################################################################################################################################
##Counts
#########################################################################################################################################
cortical_structures_thickness_count <- cortical_structures_snps2[!is.na(cortical_structures_snps2$sex) & !is.na(cortical_structures_snps2$age_imaging) & !is.na(cortical_structures_snps2$PC1) & !is.na(cortical_structures_snps2$mean_thickness_standardised)]
cortical_structures_area_count <- cortical_structures_snps2[!is.na(cortical_structures_snps2$sex) & !is.na(cortical_structures_snps2$age_imaging) & !is.na(cortical_structures_snps2$PC1) & !is.na(cortical_structures_snps2$total_surface_area_standardised)]
subcortical_structures_snps_count <- subcortical_structures_snps2[!is.na(subcortical_structures_snps2$sex) & !is.na(subcortical_structures_snps2$age_imaging) & !is.na(subcortical_structures_snps2$PC1) & !is.na(subcortical_structures_snps2$EstimatedTotalIntraCranialVol_standardised)]


descriptive_stats_cortical_thickness <- pairwiseCount(cortical_structures_thickness_count[,c(2:154)],cortical_structures_thickness_count[,c(163,165,167,169,171,173,175,177,179,181,183,185,187,189,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230)])
descriptive_stats_cortical_thickness <- t(descriptive_stats_cortical_thickness)
descriptive_stats_cortical_thickness_reshaped <-setNames(melt(descriptive_stats_cortical_thickness),c('Outcome','SNP','N'))

descriptive_stats_cortical_area <- pairwiseCount(cortical_structures_area_count[,c(2:154)],cortical_structures_area_count[,c(162,164,166,168,170,172,174,176,178,180,182,184,186,188,191,193,195,197,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227,229,231)])
descriptive_stats_cortical_area <- t(descriptive_stats_cortical_area)
descriptive_stats_cortical_area_reshaped <-setNames(melt(descriptive_stats_cortical_area),c('Outcome','SNP','N'))


descriptive_stats_subcortical <- pairwiseCount(subcortical_structures_snps_count[,c(2:154)],subcortical_structures_snps_count[,c(162:172)])
descriptive_stats_subcortical <- t(descriptive_stats_subcortical)
descriptive_stats_subcortical_reshaped <-setNames(melt(descriptive_stats_subcortical),c('Outcome','SNP','N'))

write.table(descriptive_stats_cortical_thickness_reshaped,"~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_cortical_thickness_imagen_formatted.csv",row.names=FALSE, quote=FALSE,sep=",")
write.table(descriptive_stats_cortical_area_reshaped,"~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_cortical_area_imagen_formatted.csv",row.names=FALSE, quote=FALSE,sep=",")
write.table(descriptive_stats_subcortical_reshaped,"~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_subcortical_imagen_formatted.csv",row.names=FALSE, quote=FALSE,sep=",")

#############################################################################################################################################
##RUN REGRESSION MODELS
############################################################################################################################################

#############################################################################################################################################
#CORTICAL THICKNESS
##############################################################################################################################################


setwd("~/PRS_EDU_COG/Data/IMAGEN/Genetic/Intelligence/GWAS")
grep("thickness_standardised", colnames(cortical_structures_snps2))
ivs_vec <- names(cortical_structures_snps2)[c(2:154)]
dvs_vec <- names(cortical_structures_snps2)[c(163,165,167,169,171,173,175,177,179,181,183,185,187,189,190,192,194,196,198,200,202,204,206,208,210,212,214,216,218,220,222,224,226,228,230)]
covs_vec <- paste("age_imaging","+","mean_thickness_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")


ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=cortical_structures_snps2)
})

############################################################################################################################################
#Cleaning results, extracting coefs, etc
############################################################################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)
############################################################################################################################################
## merge with sample sizes and subset results just for MR
##############################################################################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=F) 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
mr_dataset_freqs <-  merge(mr_dataset,IMAGEN_frequencies,by="SNP")
mr_dataset_freqs$eaf.outcome <- ifelse(mr_dataset_freqs$A1.x==mr_dataset_freqs$A1.y,mr_dataset_freqs$MAF,1-mr_dataset_freqs$MAF)
mr_dataset_freqs <- mr_dataset_freqs[,c(1,7,8,2,3,4,6,14)]
colnames(mr_dataset_freqs)[2:3] <- c("A1","A2")
mr_dataset_freqs<- mr_dataset_freqs[order(mr_dataset_freqs$Outcome),]
sample <- read.csv("~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_cortical_thickness_imagen_formatted.csv")
sample$SNP <- sub("_[^_]+$", "",sample$SNP)
mr_dataset_freqs <- merge(mr_dataset_freqs,sample,by=c("SNP","Outcome"))
by(mr_dataset_freqs, mr_dataset_freqs$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_T1.txt"),quote = F,row.names = F))

#############################################################################################################################################
#CORTICAL AREA
##############################################################################################################################################


grep("area_standardised", colnames(cortical_structures_snps2))
ivs_vec <- names(cortical_structures_snps2)[c(2:154)]
dvs_vec <- names(cortical_structures_snps2)[c(162,164,166,168,170,172,174,176,178,180,182,184,186,188,191,193,195,197,199,201,203,205,207,209,211,213,215,217,219,221,223,225,227,229,231)]
covs_vec <- paste("age_imaging","+","total_surface_area_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")


ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=cortical_structures_snps2)
})

############################################################################################################################################
#Cleaning results, extracting coefs, etc
############################################################################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)
############################################################################################################################################
## merge with sample sizes and subset results just for MR
##############################################################################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=F) 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
mr_dataset_freqs <-  merge(mr_dataset,IMAGEN_frequencies,by="SNP")
mr_dataset_freqs$eaf.outcome <- ifelse(mr_dataset_freqs$A1.x==mr_dataset_freqs$A1.y,mr_dataset_freqs$MAF,1-mr_dataset_freqs$MAF)
mr_dataset_freqs <- mr_dataset_freqs[,c(1,7,8,2,3,4,6,14)]
colnames(mr_dataset_freqs)[2:3] <- c("A1","A2")
mr_dataset_freqs<- mr_dataset_freqs[order(mr_dataset_freqs$Outcome),]
sample <- read.csv("~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_cortical_area_imagen_formatted.csv")
sample$SNP <- sub("_[^_]+$", "",sample$SNP)
mr_dataset_freqs <- merge(mr_dataset_freqs,sample,by=c("SNP","Outcome"))
by(mr_dataset_freqs, mr_dataset_freqs$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_T1.txt"),quote = F,row.names = F))



#############################################################################################################################################
#SUBCORTICAL VOLUMES
##############################################################################################################################################


ivs_vec <- names(subcortical_structures_snps2)[c(2:154)]
dvs_vec <- names(subcortical_structures_snps2)[c(162:172)]
covs_vec <- paste("age_imaging","+","EstimatedTotalIntraCranialVol_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")


ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=subcortical_structures_snps2)
})

############################################################################################################################################
#Cleaning results, extracting coefs, etc
############################################################################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)
############################################################################################################################################
## merge with sample sizes and subset results just for MR
##############################################################################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=F) 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
mr_dataset_freqs <-  merge(mr_dataset,IMAGEN_frequencies,by="SNP")
mr_dataset_freqs$eaf.outcome <- ifelse(mr_dataset_freqs$A1.x==mr_dataset_freqs$A1.y,mr_dataset_freqs$MAF,1-mr_dataset_freqs$MAF)
mr_dataset_freqs <- mr_dataset_freqs[,c(1,7,8,2,3,4,6,14)]
colnames(mr_dataset_freqs)[2:3] <- c("A1","A2")
mr_dataset_freqs<- mr_dataset_freqs[order(mr_dataset_freqs$Outcome),]
sample <- read.csv("~/PRS_EDU_COG/Data/IMAGEN/Intelligence/Descriptive_stats/descriptive_stats_subcortical_imagen_formatted.csv")
sample$SNP <- sub("_[^_]+$", "",sample$SNP)
mr_dataset_freqs <- merge(mr_dataset_freqs,sample,by=c("SNP","Outcome"))
by(mr_dataset_freqs, mr_dataset_freqs$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_T1.txt"),quote = F,row.names = F))


