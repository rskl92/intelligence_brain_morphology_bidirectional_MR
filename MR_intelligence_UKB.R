freqs <- fread("~/Data/UK_Biobank/GeneticData/ukbb_intelligence_freqs.frq")
colnames(freqs)[3:5] <- c("A1.freq","A2.freq","freq")
##read in exposure
exposure <- read_exposure_data("~/Data/Discovery_samples/Intelligence/cognitive_ability_gwas_snps_full.txt", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", pval_col = "P")
exposure<- merge(exposure,freqs,by="SNP")
exposure$eaf.exposure <- ifelse(exposure$effect_allele.exposure==exposure$A1.freq,exposure$freq,1-exposure$freq)

setwd("~~/Data/UK_Biobank/GWAS/Intelligence/")
filelist<- list.files()
l<- list()
for(i in filelist)
{
  l[[i]] <- read_outcome_data(filename=i,snp_col="SNP",phenotype_col="Outcome",beta_col ="BETA",se_col="SE",pval_col="P",samplesize_col="N")
l[[i]]$Outcome <- i
}
outcome_dat <- bind_rows(l)
outcome_dat$Outcome<- gsub(pattern = "\\.txt$", "", outcome_dat$Outcome)
outcome_dat$A1 <- sub("^[^_]*_", "",outcome_dat$SNP)
outcome_dat$SNP <- sub("_[^_]+$","",outcome_dat$SNP)
outcome_dat$A1 <-toupper(outcome_dat$A1)
merged_outcome_exp <- merge(outcome_dat,exposure,by="SNP")
merged_outcome_exp$A2 <- ifelse(merged_outcome_exp$A1==merged_outcome_exp$A1.freq,merged_outcome_exp$A2.freq,merged_outcome_exp$A1.freq)
merged_outcome_exp$eaf.outcome <- ifelse(merged_outcome_exp$A1==merged_outcome_exp$A1.freq,merged_outcome_exp$freq,1-merged_outcome_exp$freq)
merged_outcome_exp <- merged_outcome_exp[,c(1:5,14:15,32)]
write.table(merged_outcome_exp,"~/Data/UK_Biobank/UK_Biobank_GWAs_intelligence_brain.csv",sep=",",row.names=FALSE,quote=FALSE)



outcome_dat <- read_outcome_data("~/Data/UK_Biobank/UK_Biobank_GWAs_intelligence_brain.csv",sep=",",phenotype_col = "Outcome",beta_col ="beta.outcome",se_col="se.outcome", eaf_col="eaf.outcome", effect_allele_col = "A1", other_allele_col = "A2",samplesize_col = "samplesize.outcome")
dat <- harmonise_data(exposure,outcome_dat)


#Running the MR
mr_results <- mr(dat)
mr_results <- mr_results[!mr_results$method=="simple mode",]

#Set results directory
setwd("~/Results/UK_Biobank/Intelligence")

#Putting MR results into a HTML report
mr_report(dat, output_path =paste("Alzheimer's disease","dat$brainStructure"), author="Analyst", study = paste("Alzheimer's disease","dat$brainStructure",sep=""))

#Save results
write.csv(mr_results,"mr_results_cognitive_to_all_brainstructures.txt",quote=FALSE,row.names = FALSE)



#Prepare columns for steiger filtering
dat<- harmonise_data(exposure,outcome_dat)
dat <- dat[dat$mr_keep=="TRUE",]
dat$units.exposure <- "SD"
dat$units.outcome <- "SD"
dat$samplesize.exposure=248482


#Perform steiger
dat2<- steiger_filtering(dat)
directionality <- directionality_test(dat2)
write.table(dat2,"steiger_filtering_intelligence_to_all_brain_structures.csv", row.names=F,quote=F)
write.table(directionality,"directionality_intelligence_to_all_brain_structures.csv", row.names=F,quote=F)


#pleiotropy test
pleiotropy <- mr_pleiotropy_test(dat)
pleiotropy <- pleiotropy[!duplicated(pleiotropy),]
write.csv(pleiotropy,"pleiotropy_mr.csv",row.names=F,quote=F)

#heterogeneity test
heterogeneity <- mr_heterogeneity(dat)
heterogeneity <- heterogeneity[heterogeneity$method=="Inverse variance weighted",]
write.csv(heterogeneity,"heterogeneity_mr.csv",row.names=FALSE,quote=FALSE)

#Save dataframe
write.csv(dat2, file="full_dat.csv", row.names = FALSE)


