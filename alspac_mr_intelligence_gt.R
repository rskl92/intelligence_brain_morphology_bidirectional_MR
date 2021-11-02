install.packages("devtools")
install.packages("digest")
install.packages("rlang")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
devtools::install_github("MRCIEU/xlsx")
install.packages("xlsx")
install.packages("psych")
install.packages("plyr")
install.packages("ggplot2")
install.packages("microbenchmark")

## This calls the packages you just installed so that you can use them ##
library(rlang)
library(digest)
library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(xlsx)
library(psych)
library(ggplot2)
library(plyr)
library(dplyr)
library(microbenchmark)
##############################################################################################################################
########################################## MR  ################################################################
##############################################################################################################################

##############################################################################################################################
#Setting your working directory where everything will be called from and also saved
##############################################################################################################################

rm(list=ls())

ukb_intelligence_freqs<- fread("~/PRS_EDU_COG/Data/UKBiobank/Genetic/ukb_intelligence_instruments_and_proxies_frqs.afreq")
exposure_dat<- read_exposure_data("~/PRS_EDU_COG/Data/Discovery_Samples/Intelligence/intelligence_instruments_proxies020621_exposure_MR.txt",snp_col = "SNP", beta_col = "Beta", se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P")
exposure_dat <- merge(exposure_dat,ukb_intelligence_freqs,by.x="SNP",by.y="ID")
exposure_dat$eaf.exposure<- ifelse(exposure_dat$effect_allele.exposure==exposure_dat$ALT,exposure_dat$ALT_FREQS,1-exposure_dat$ALT_FREQS)
exposure_dat <- exposure_dat[,c(1:12)]

##############################################################################################################################
#naming exposure 
##############################################################################################################################

exposure_dat$exposure<- "Intelligence"


##############################################################################################################################
##set to directory with ALSPAC GWAS
##############################################################################################################################

setwd("~/PRS_EDU_COG/Data/ALSPAC/Genetic/Intelligence/GWAS")

##############################################################################################################################
#Read in files containing effects of intelligence SNPs on brain structures in ALSPAC
##############################################################################################################################

myfiles<- list.files(pattern=".txt")
l <- list()
for (i in myfiles) 
{
  l[[i]] <- read_outcome_data(i,snp_col="SNP",sep=" ",phenotype_col="brainStructure",eaf_col = "EAF",effect_allele_col="A1",other_allele_col="A2",beta_col ="BETA",se_col="SE",pval_col="P",samplesize_col="N")
  l[[i]]$outcome <- i
}
outcome_dat <- bind_rows(l)

##############################################################################################################################
#harmonising exposure and outcome datasets (making sure alleles aligned)
##############################################################################################################################

dat <- harmonise_data(exposure_dat, outcome_dat,action=1) 

##############################################################################################################################
#plot eaf of exposure and outcome to see if they are on the same strand
##############################################################################################################################

plot_eafs <- plot(dat$eaf.outcome,dat$eaf.exposure)

##############################################################################################################################
#Running the MR
##############################################################################################################################

mr_results <- mr_results <- mr(dat)


##############################################################################################################################
#Remove results using simple mode
##############################################################################################################################


mr_results <- mr_results[!mr_results$method=="simple mode",]

##############################################################################################################################
#displayng MR results 
##############################################################################################################################

mr_results

##############################################################################################################################
#Set results directory
##############################################################################################################################

setwd("~/PRS_EDU_COG/Results/Intelligence/ALSPAC")

##############################################################################################################################
#Putting MR results into a HTML report
##############################################################################################################################

write.csv(mr_results,"mr_results_intelligence_alspac_brain.txt",quote=F,row.names = F)


##############################################################################################################################
#Putting MR results into a HTML report
##############################################################################################################################

mr_report(dat, output_path = "cortical_brain_alspac.txt", author="Analyst", study = paste("intelligence","-brain",sep=""))

##############################################################################################################################
##Perform steiger filtering
##############################################################################################################################


dat$units.exposure <- "SD"
dat$units.outcome <- "SD"
dat$samplesize.exposure=248482


##############################################################################################################################
# Calculate F statistics
##############################################################################################################################

dat$F   = dat$BXG^2/dat$seBetaXG^2

##############################################################################################################################
#Sensitivity analyses
##############################################################################################################################

directionality <- directionality_test(dat)
pleiotropy <- mr_pleiotropy_test(dat)
heterogeneity <- mr_heterogeneity(dat)
steiger <- steiger_filtering(dat)

##############################################################################################################################
#Store results in directory
##############################################################################################################################

write.csv(heterogeneity,"heterogeneity_alspac_intelligence_all_brain_structures.csv",row.names=F,quote=F)
write.csv(pleiotropy,"pleiotropy_alspac_intelligence_all_brain_structures.csv",row.names=F,quote=F)
write.csv(directionality,"directionality_alspac_intelligence_all_brain_structures.csv",row.names=F,quote=F)
write.csv(steiger, file="full_dat_intelligence_brain.csv", row.names = FALSE,quote=FALSE)