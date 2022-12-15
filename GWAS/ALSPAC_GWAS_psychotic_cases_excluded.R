install.packages(repurrrsive)
install.packages(purrr)
install.packages(stringr)
install.packages(broom)
install.packages(broom)
install.packages(data.table)
install.packages(readr)
install.packages(tidyverse)
install.packages(readxl)
install.packages(psych)
install.packages(data.table)
install.packages(naniar)


library(haven)
library(reshape2)
library(dplyr)
library(repurrrsive)
library(purrr)
library(stringr)
library(broom)
library(broom)
library(data.table)
library(readr)
library(tidyverse)
library(readxl)
library(psych)
library(data.table)


##WRITE FUNCTIONS TO BE USED IN SCRIPTS

colClean <- function(x){ colnames(x) <- gsub("thic.+$", "thickness", colnames(x))
colnames(x) <- gsub("age_months", "age", colnames(x))
colnames(x) <- gsub("lthickness","l_mean_thickness",colnames(x))
colnames(x) <- gsub("rthickness","r_mean_thickness",colnames(x))
colnames(x) <- gsub("surf.+$", "area", colnames(x))
colnames(x) <- gsub("larea","l_total_surface_area",colnames(x))
colnames(x) <- gsub("rarea","r_total_surface_area",colnames(x))
colnames(x) <- gsub("laccumb","l_Accumbens",colnames(x))
colnames(x) <- gsub("raccumb","r_Accumbens",colnames(x))
colnames(x) <- gsub("lpal","l_Pallidum",colnames(x))
colnames(x) <- gsub("rpal","r_Pallidum",colnames(x))
colnames(x) <- gsub("lthal","l_Thalamus",colnames(x))
colnames(x) <- gsub("rthal","r_Thalamus",colnames(x))
colnames(x) <- gsub("lput","l_Putamen",colnames(x))
colnames(x) <- gsub("rput","r_Putamen",colnames(x))
colnames(x) <- gsub("rhippo","r_Hippocampus",colnames(x))
colnames(x) <- gsub("lhippo","l_Hippocampus",colnames(x))
colnames(x) <- gsub("lcaud","l_Caudate",colnames(x))
colnames(x) <- gsub("rcaud","r_Caudate",colnames(x))
colnames(x) <- gsub("rlatvent","r_LateralVentricle",colnames(x))
colnames(x) <- gsub("llatvent","l_LateralVentricle",colnames(x))
colnames(x) <- gsub("icv","EstimatedTotalIntraCranialVol",colnames(x))
colnames(x) <- gsub("lamyg","l_Amygdala",colnames(x))
colnames(x) <- gsub("ramyg","r_Amygdala",colnames(x)); x }



##Read in SNPs 
intelligence_snps <- fread("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Intelligence/SNPs/alspac_intelligence_hill.raw")
alspac_freqs_intelligence <- fread("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Intelligence/SNPs/ALSPAC_freqs_intelligence_hill.frq")
intelligence_snps$aln <- as.numeric(str_extract(intelligence_snps$FID, "[0-9]+"))
intelligence_snps$qlet <- (str_extract(intelligence_snps$FID, "[aA-zZ]+"))

##Read in genetic principal components
PCs <- read_table2("~/ALSPAC_brain/PCs/data.eigenvec", col_names = FALSE)
PCs <- PCs[,2:7]
names(PCs)[-1]<- paste0('PC',1:(ncol(PCs)-1))
PCs$aln<- as.numeric(str_extract(PCs$X2, "[0-9]+"))
PCs$qlet <- (str_extract(PCs$X2, "[aA-zZ]+"))



##Merge genetic and phenotypic datasets
merge_genetic_intelligence <- merge(intelligence_snps,PCs,by=c("aln","qlet"))

##Import MRI data
ALSPAC_MRI_extraction_data_alnqlet <- read_dta("~/ALSPAC_brain/ALSPAC-MRI_extraction_data_alnqlet.dta")
ALSPAC_MRI_extraction_data_alnqlet <- colClean(ALSPAC_MRI_extraction_data_alnqlet)

##3 B individuals, 2 of which do not have an A- remove 1 from pair row 50
which(ALSPAC_MRI_extraction_data_alnqlet$aln==31365 & ALSPAC_MRI_extraction_data_alnqlet$qlet=="B")
ALSPAC_MRI_extraction_data_alnqlet<- ALSPAC_MRI_extraction_data_alnqlet[-50,]
ALSPAC_MRI_extraction_data_alnqlet[ALSPAC_MRI_extraction_data_alnqlet$aln==33082,] <- NA

##Keep controls of the ALSPAC psychosis substudy
ALSPAC_MRI_extraction_data_alnqlet<- ALSPAC_MRI_extraction_data_alnqlet[ALSPAC_MRI_extraction_data_alnqlet$PE_status==0 | is.na(ALSPAC_MRI_extraction_data_alnqlet$PE_status),] 


##Recoding values to missing
ALSPAC_MRI_extraction_data_alnqlet[ALSPAC_MRI_extraction_data_alnqlet==-1] <- NA
ALSPAC_MRI_extraction_data_alnqlet[ALSPAC_MRI_extraction_data_alnqlet==-9999] <- NA
ALSPAC_MRI_extraction_data_alnqlet[ALSPAC_MRI_extraction_data_alnqlet==0] <- NA
ALSPAC_MRI_extraction_data_alnqlet$sex <- as.factor(ALSPAC_MRI_extraction_data_alnqlet$sex)

##subset dataset into cortical and subcortical for separate qc filtering (i.e. for subcortical structures data keep the those imaging measures and the variable denoting subcortical QC rating)
subcortical_structures <- ALSPAC_MRI_extraction_data_alnqlet[,c(1:11,84,153:168,170)]

##Keep cortical structures by generating a new dataset without subcortical structures
cortical_structures <- ALSPAC_MRI_extraction_data_alnqlet[,-c(84,153:168,170)]

##FILTER BASED ON QC (I.E. filter subcortical structures based on subcortical qc rating and cortical structures based on cortical QC rating)
cortical_structures[which(cortical_structures$corticalqc==2),c(12:79,84:151)]<- NA
subcortical_structures[which(subcortical_structures$subcorticalqc==2),c(13:28)] <- NA


##select columns to put into subcortical dataset and the cortical dataset(i.e. subcortical volumes and cortical thickness and surface area)
subcortical_structures2 <- subcortical_structures[,c(12:28)]
cortical_structures2 <- cortical_structures[,c(12:151)]

#################################################################################################
##calculate mean for cortical structures
#################################################################################################

out_cort<- t(cortical_structures2) %>%
  data.frame() %>%
  group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
  summarise_all(mean) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t()


#################################################################################################
##add ID to column
#################################################################################################

out2_cort <- cbind(cortical_structures[,c(1:11)], out_cort)
save(out2_cort,file="cortical_structures_unstandardised_clean.Rdata")

#################################################################################################
##standardise mean measures
#################################################################################################

data.standardised_cortical <- as.data.frame(sapply(out2_cort[,c(12:81)], scale))
names(data.standardised_cortical) <- sapply(names(data.standardised_cortical), function(x) {
  y <- paste0(x, "_standardised")
}
)

#################################################################################################
##bind standardised cortical measures with main cortical dataset
#################################################################################################

cortical_structures <- cbind(out2_cort[1:11],data.standardised_cortical)



#################################################################################################
##calculate mean for subcortical structures
#################################################################################################

out_subcort<-t(subcortical_structures2) %>%
  data.frame() %>%
  group_by(., id = gsub('^[rl]_', '', rownames(.))) %>% ##replace beginning that starts with r or l with nothing and group by means group them if they have the same name
  summarise_all(mean) %>%
  data.frame() %>%
  column_to_rownames(var = 'id') %>%
  t()


#################################################################################################
##add IDs (aln and qlet)
#################################################################################################

out2_subcort <- cbind(subcortical_structures[,c(1:11)], out_subcort)
save(out2_subcort,file="subcortical_structures_unstandardised_clean.Rdata")


#################################################################################################
##standardise mean measures
#################################################################################################

data.standardised_subcortical <- as.data.frame(sapply(out2_subcort[,12:20], scale))
names(data.standardised_subcortical) <- sapply(names(data.standardised_subcortical), function(x) {
  y <- paste0(x, "_standardised")
}
)



#################################################################################################
##bind standardised measures with main cortical dataset
#################################################################################################

subcortical_structures <- cbind(out2_subcort[1:11],data.standardised_subcortical)


##merge subcortical and genetic data
intelligence_genetic_subcortical_structures<- merge(merge_genetic_intelligence,subcortical_structures,by=c("aln","qlet"))
intelligence_genetic_cortical_structures <- merge(merge_genetic_intelligence,cortical_structures,by=c("aln","qlet"))

##separate cortical area and thickness for intelligence SNPs into separate datasets
paste(grep("thick",colnames(intelligence_genetic_cortical_structures)),collapse=",")
paste(grep("area",colnames(intelligence_genetic_cortical_structures)),collapse=",")
intelligence_genetic_cortical_area <- intelligence_genetic_cortical_structures[,c(1:177,179,181,183,185,187,189,191,193,195,197,199,201,203,206,208,210,212,214,216,218,220,222,224,226,228,230,232,234,236,238,240,242,244,245)]
intelligence_genetic_cortical_thickness <- intelligence_genetic_cortical_structures[,c(1:176,178,180,182,184,186,188,190,192,194,196,198,200,202,204,205,207,209,211,213,215,217,219,221,223,225,227,229,231,233,235,237,239,241,243,246)]


##filter dataset to include complete data to be used in regression models
intelligence_genetic_cortical_thickness <- intelligence_genetic_cortical_thickness[!is.na(intelligence_genetic_cortical_thickness$sex) & !is.na(intelligence_genetic_cortical_thickness$mean_thickness_standardised) & !is.na(intelligence_genetic_cortical_thickness$PC1) & !is.na(intelligence_genetic_cortical_area$age)]
intelligence_genetic_cortical_area <- intelligence_genetic_cortical_area[!is.na(intelligence_genetic_cortical_area$sex) & !is.na(intelligence_genetic_cortical_area$total_surface_area_standardised) & !is.na(intelligence_genetic_cortical_area$age)  & !is.na(intelligence_genetic_cortical_area$PC1)]
intelligence_genetic_subcortical_structures <- intelligence_genetic_subcortical_structures[!is.na(intelligence_genetic_subcortical_structures$sex) & !is.na(intelligence_genetic_subcortical_structures$EstimatedTotalIntraCranialVol_standardised)  & !is.na(intelligence_genetic_subcortical_structures$PC1)]

intelligence_descriptive_stats_cortical_area <- pairwiseCount(intelligence_genetic_cortical_area[,c(9:161)],intelligence_genetic_cortical_area[,c(177:211)])
intelligence_descriptive_stats_cortical_area <- t(intelligence_descriptive_stats_cortical_area)
intelligence_descriptive_stats_cortical_reshaped_area <-setNames(melt(intelligence_descriptive_stats_cortical_area),c('Outcome','SNP','N'))
write.table(intelligence_descriptive_stats_cortical_reshaped_area,"~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_cortical_area_alspac_sample.csv",row.names=FALSE, quote=FALSE,sep=",")

intelligence_descriptive_stats_cortical_thickness <- pairwiseCount(intelligence_genetic_cortical_thickness[,c(9:161)],intelligence_genetic_cortical_thickness[,c(177:211)])
intelligence_descriptive_stats_cortical_thickness <- t(intelligence_descriptive_stats_cortical_thickness)
intelligence_descriptive_stats_cortical_reshaped_thickness <-setNames(melt(intelligence_descriptive_stats_cortical_thickness),c('Outcome','SNP','N'))
write.table(intelligence_descriptive_stats_cortical_reshaped_thickness,"~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_cortical_thickness_alspac_sample.csv",row.names=FALSE, quote=FALSE,sep=",")

intelligence_descriptive_stats_subcortical_structures <- pairwiseCount(intelligence_genetic_subcortical_structures[,c(9:161)],intelligence_genetic_subcortical_structures[,c(177:185)])
intelligence_descriptive_stats_subcortical_structures <- t(intelligence_descriptive_stats_subcortical_structures)
intelligence_descriptive_stats_subcortical_structures_reshaped <-setNames(melt(intelligence_descriptive_stats_subcortical_structures),c('Outcome','SNP','N'))
write.table(intelligence_descriptive_stats_subcortical_structures_reshaped,"~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_subcortical_alspac_sample.csv",row.names=FALSE, quote=FALSE,sep=",")

##remove (/ from rsids) or else regression models will show an error because of specia characters
names(intelligence_genetic_cortical_area) <- gsub("\\(..)", "",names(intelligence_genetic_cortical_area))
names(intelligence_genetic_cortical_thickness) <- gsub("\\(..)", "",names(intelligence_genetic_cortical_thickness))
names(intelligence_genetic_subcortical_structures) <- gsub("\\(..)", "", names(intelligence_genetic_subcortical_structures))

# Save dataset (only saving intelligence as it is a smaller file. The samples should be the same)
save(intelligence_genetic_cortical_area,file="intelligence_genetic_cortical_area.Rdata")
save(intelligence_genetic_cortical_thickness,file="intelligence_genetic_cortical_thickness.Rdata")
save(intelligence_genetic_subcortical_structures,file="intelligence_genetic_cortical_subcortical_structures.Rdata")





###################################################################################################
#SET WORKING DIRECTORY TO OUTPUT THE BRAIN GWAS FOR THE INTELLIGENCE SNPS
###################################################################################################

setwd("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Intelligence/GWAS")


###################################################################################################
#GWAS OF INTELLIGENCE SNPS ON REGIONAL AND GLOBAL THICKNESS
###################################################################################################


ivs_vec <- names(intelligence_genetic_cortical_thickness)[9:161]
dvs_vec <- names(intelligence_genetic_cortical_thickness)[c(177:211)]
covs_vec <- paste("age","+","mean_thickness_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")

ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=intelligence_genetic_cortical_thickness)
})

# Creating / combining results
tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)


# subset results just for MR


x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=FALSE); 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
intelligence_alleles <- read.table("~/PRS_EDU_COG/Data/Discovery_Samples/Intelligence/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
sample_size <- read.csv("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_cortical_thickness_alspac_sample.csv")
sample_size$SNP<- sub("_[^_]+$", "",sample_size$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
colnames(mr_dataset) <- c("SNP","Outcome","BETA","SE","T-statistic","P","A1","A2")
mr_dataset <- merge(mr_dataset,sample_size,by.x=c("SNP","Outcome"),by.y=c("SNP","Outcome"))
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
mr_dataset <- merge(mr_dataset,alspac_freqs_intelligence,by="SNP")
mr_dataset$EAF <- ifelse(mr_dataset$A1.x==mr_dataset$A1.y,mr_dataset$MAF,1-mr_dataset$MAF)
mr_dataset <- mr_dataset[,c(1,2,7:8,3:6,9,10,15)]
colnames(mr_dataset) <- c("SNP","Outcome","A1","A2","BETA","SE","T-statistic","P","N","CHR","EAF")
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_ALSPAC_without_pe_cases.txt"),quote = FALSE,row.names = FALSE))



#####################################################################################################
#GWAS OF INTELLIGENCE SNPs ON REGIONAL AND GLOBAL SURFACE AREA
###################################################################################################

ivs_vec <- names(intelligence_genetic_cortical_area)[9:161]
dvs_vec <- names(intelligence_genetic_cortical_area)[c(177:211)]
covs_vec <- paste("age","+","total_surface_area_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")

###################################################################################################
##CREATE FORMULAS FOR LINEAR REFRESSION MODELS
###################################################################################################
ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=intelligence_genetic_cortical_area)
})

###################################################################################################
#COMBINED THE RESULTS FROM THE REGRESSION MODELS
###################################################################################################

tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)

###################################################################################################
# subset results just for MR
###################################################################################################

x <- lapply(tidy_results,"[",2,,drop=FALSE) ##2ND COLUMN REFERS TO THE B-COEFFICIENT
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names=FALSE); 
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
intelligence_alleles <- read.table("~/PRS_EDU_COG/Data/Discovery_Samples/Intelligence/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
sample_size <- read.csv("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_cortical_area_alspac_sample.csv")
sample_size$SNP<- sub("_[^_]+$", "",sample_size$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
colnames(mr_dataset) <- c("SNP","Outcome","BETA","SE","T-statistic","P","A1","A2")
mr_dataset <- merge(mr_dataset,sample_size,by.x=c("SNP","Outcome"),by.y=c("SNP","Outcome"))
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
mr_dataset <- merge(mr_dataset,alspac_freqs_intelligence,by="SNP")
mr_dataset$EAF <- ifelse(mr_dataset$A1.x==mr_dataset$A1.y,mr_dataset$MAF,1-mr_dataset$MAF)
mr_dataset <- mr_dataset[,c(1,2,7:8,3:6,9,10,15)]
colnames(mr_dataset) <- c("SNP","Outcome","A1","A2","BETA","SE","T-statistic","P","N","CHR","EAF")
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_ALSPAC_without_pe_cases.txt"),quote = FALSE,row.names = FALSE))



####################################################################################################
##GWAS OF INTELLIGENCE SNPs ON SUBCORTICAL VOLUMES
####################################################################################################

ivs_vec <- names(intelligence_genetic_subcortical_structures)[9:161]
dvs_vec <- names(intelligence_genetic_subcortical_structures)[177:185]
covs_vec <- paste("age","+","EstimatedTotalIntraCranialVol_standardised","+","sex","+","PC1","+","PC2","+","PC3","+","PC4","+","PC5")

ivs <- paste0(" ~ ", ivs_vec,"+",covs_vec)
dvs_ivs <- unlist(lapply(ivs, function(x) paste0(dvs_vec, x))) ##given a list x, unlist simplifies to produce a vector
formulas <- lapply(dvs_ivs, formula)

lm_results <- lapply(formulas, function(x) {
  lm(x, data=intelligence_genetic_subcortical_structures)
})


# Creating / combining results
tidy_results <- lapply(lm_results, broom::tidy)
dv_list <- as.list(stringi::stri_extract_first_words(dvs_ivs))
tidy_results <- Map(cbind, dv_list, tidy_results)
# subset results just for MR

x <- lapply(tidy_results,"[",2,,drop=FALSE)
rs_identifiers <- Map(cbind,dv_list,x)
rs <- Map(as.data.frame,rs_identifiers)
dfrData_rs <- rbindlist(rs,use.names = FALSE)
dfrData_rs$`dots[[1L]][[1L]]`<-NULL
colnames(dfrData_rs) <- c("Outcome", "SNP", "BETA", "SE", "Statistic", "P")
combined_results <- dfrData_rs
combined_results$SNP <- sub("_[^_]+$", "",combined_results$SNP)
intelligence_alleles <- read.table("~/PRS_EDU_COG/Data/Discovery_Samples/Intelligence/intelligence_instruments_proxies020621_exposure_gwas_A1_A2.txt",he=TRUE)
sample_size <- read.csv("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Descriptive_stats/Intelligence/intelligence_descriptive_stats_subcortical_alspac_sample.csv")
sample_size$SNP<- sub("_[^_]+$", "",sample_size$SNP)
mr_dataset <- merge(combined_results,intelligence_alleles,by="SNP")
colnames(mr_dataset) <- c("SNP","Outcome","BETA","SE","T-statistic","P","A1","A2")
mr_dataset <- merge(mr_dataset,sample_size,by.x=c("SNP","Outcome"),by.y=c("SNP","Outcome"))
mr_dataset<- mr_dataset[order(mr_dataset$Outcome),]
mr_dataset <- merge(mr_dataset,alspac_freqs_intelligence,by="SNP")
mr_dataset$EAF <- ifelse(mr_dataset$A1.x==mr_dataset$A1.y,mr_dataset$MAF,1-mr_dataset$MAF)
mr_dataset <- mr_dataset[,c(1,2,7:8,3:6,9,10,15)]
colnames(mr_dataset) <- c("SNP","Outcome","A1","A2","BETA","SE","T-statistic","P","N","CHR","EAF")
by(mr_dataset, mr_dataset$Outcome, FUN=function(i) write.table(i, paste0(i$Outcome[1], "_ALSPAC_without_pe_cases.txt"),quote = FALSE,row.names = FALSE))


