# ml R/4.2.1-foss-2020b
#################################################5e-8
library(TwoSampleMR)
library(dplyr)
'%ni%' <- Negate('%in%')
library(MendelianRandomization)
library(MRPRESSO)

ao <- available_outcomes()

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_met_MS_06_22_22/"

####### met-a
data_met_a=ao
data_met_a$new_name_1=sapply(strsplit(data_met_a$id, split= "-"),"[",1)
data_met_a$new_name_2=sapply(strsplit(data_met_a$id, split= "-"),"[",2)
data_met_a=data_met_a[data_met_a$new_name_1=="met",]
data_met_a=data_met_a[data_met_a$new_name_2=="a",]

####### met-c
data_met_c=ao
data_met_c$new_name_1=sapply(strsplit(data_met_c$id, split= "-"),"[",1)
data_met_c$new_name_2=sapply(strsplit(data_met_c$id, split= "-"),"[",2)
data_met_c=data_met_c[data_met_c$new_name_1=="met",]
data_met_c=data_met_c[data_met_c$new_name_2=="c",]

####### met-d
data_met_d=ao
data_met_d$new_name_1=sapply(strsplit(data_met_d$id, split= "-"),"[",1)
data_met_d$new_name_2=sapply(strsplit(data_met_d$id, split= "-"),"[",2)
data_met_d=data_met_d[data_met_d$new_name_1=="met",]
data_met_d=data_met_d[data_met_d$new_name_2=="d",]

#b_Wald, se_Wald, pval_Wald,
final_res=data.frame("id.exposure","exposure","id.outcome","outcome",
                     "b_IVW_MRE", "se_IVW_MRE", "pval_IVW_MRE", "b_IVW_FE", "se_IVW_FE", "pval_IVW_FE",
                     "b_Egger", "se_Egger", "pval_Egger",
                     "Egger_intercept", "pval_intercept", "Het_IVW_pval", "Het_Egger_pval",
                     "b_W_Med", "se_W_Med", "pval_W_Med", "b_W_Mod", "se_W_Mod", "pval_W_Mod",
                     "nsnps","Total_R_Square", "F_stat", "F_stat_sim", "F_statistic",
                     "snp_r2.exposure","snp_r2.outcome","causal_direction","steiger_pval",
                     "b_PRESSO_raw", "se_PRESSO_raw", "pval_PRESSO_raw", 
                     "b_PRESSO_corrected", "se_PRESSO_corrected", "pval_PRESSO_corrected", "pval_PRESSO_Global")

write.table(final_res, file= paste(Pathway_out,"MenRan_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')
# write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
#             row.names = F, quote = FALSE, na = "-",sep='\t')


all_data_met <- c(data_met_a$id, data_met_c$id)
# all_data_met <- all_data_met[-21]
# all_data_met <- all_data_met[-20]
# all_data_met <- all_data_met[-19]
# 
# print(all_data_met)
options(scipen = 5)

for (e in all_data_met) {
  for (n in c("ieu-b-18")) {
    try(exp_dat <- extract_instruments(outcomes=e, p1=5e-8),silent = T)
    if (exists("exp_dat")==TRUE){
      if (length(exp_dat$SNP)>0) {
        exp_dat$exposure=sapply(strsplit(exp_dat$exposure,fixed = TRUE, split= " ||"),"[",1)
        outcome_dat <-extract_outcome_data(
          snps = exp_dat$SNP,
          outcomes=n)
        if (length(outcome_dat$SNP)>0) {
          outcome_dat$outcome=sapply(strsplit(outcome_dat$outcome,fixed = TRUE, split= " ||"),"[",1)
          #My_MR(exp_dat,outcome_dat)
          rm(dat)
          try(dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat), silent=TRUE)
          dat=dat[is.na(dat$beta.outcome)==FALSE ,]
          
          ############# Change
          dat=dat[dat$mr_keep==TRUE,]
          
          if (length(which(dat$mr_keep=='TRUE'))>2) {
            #MendelianRandomization R Package
            MR_IVW_MRE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "random")
            MR_IVW_FE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "fixed")
            MR_Egger = MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                 by = dat$beta.outcome, byse = dat$se.outcome))
            MR_W_Med = MendelianRandomization::mr_median(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                  by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            MR_W_Mod = mr_mbe(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                       by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            
            dat_steiger = steiger_filtering(dat)
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.exposure")] <- "r.exposure"
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.outcome")] <- "r.outcome"
            
            dat_steiger_direct = directionality_test(dat_steiger)
            
            dat$len_SNP=length(dat$SNP)
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5082560/
            try(dat$R_Square<- 2*(dat$beta.exposure^2)*dat$eaf.exposure*(1-dat$eaf.exposure), silent=TRUE)
            try(dat$Total_R_Square <- sum(dat$R_Square), silent=TRUE)
            #try(dat$F_stat_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            try(dat$Total_R_Square <- dat$Total_R_Square/(dat$SD^2), silent=TRUE)
            try(dat$F_stat <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            #try(dat$F_stat_sim_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_stat_sim <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_statistic_pre <- (dat$beta.exposure/dat$se.exposure)^2, silent=TRUE)
            try(dat$F_statistic <- sum(dat$F_statistic_pre)/(dat$len_SNP), silent=TRUE)
            try(F_stat <- data.frame(dat$id.exposure[1],dat$len_SNP[1],dat$samplesize.exposure[1],dat$Total_R_Square[1],dat$F_stat[1],dat$F_stat_sim[1],dat$F_statistic[1]), silent=TRUE)
            
            if (length(which(dat$mr_keep=='TRUE'))>3) {
              PressoObject=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                                     SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, 
                                     NbDistribution = max(1000,length(dat$SNP)/0.05),  SignifThreshold = 0.05)
              PRESSO = PressoObject$`Main MR results`
              PRESSO_raw = PRESSO[PRESSO$"MR Analysis"=='Raw',]
              PRESSO_corrected = PRESSO[PRESSO$"MR Analysis"=='Outlier-corrected',]
              PRESSO_GLOBAL = PressoObject$`MR-PRESSO results`$`Global Test`
              pval_PRESSO_Global = PRESSO_GLOBAL$Pvalue
              names(PRESSO_raw)[names(PRESSO_raw) == "Causal Estimate"] <- "b_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "Sd"] <- "se_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "P-value"] <- "pval_PRESSO_raw"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Causal Estimate"] <- "b_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Sd"] <- "se_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "P-value"] <- "pval_PRESSO_corrected"
              
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   PRESSO_raw$b_PRESSO_raw, PRESSO_raw$se_PRESSO_raw, PRESSO_raw$pval_PRESSO_raw, 
                                   PRESSO_corrected$b_PRESSO_corrected, PRESSO_corrected$se_PRESSO_corrected, 
                                   PRESSO_corrected$pval_PRESSO_corrected, pval_PRESSO_Global)
              
              write.table(final_res, file= paste(Pathway_out,"MenRan_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
              
            } else {
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   NA,NA,NA,NA,NA,NA,NA)
              
              write.table(final_res, file= paste(Pathway_out,"MenRan_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
            }
          }
        }
      }
    }
  }
}

all_data_met <- c(data_met_d$id)
# all_data_met <- all_data_met[-21]
# all_data_met <- all_data_met[-20]
# all_data_met <- all_data_met[-19]
# 
# print(all_data_met)

for (e in all_data_met) {
  for (n in c("ieu-b-18")) {
    try(exp_dat <- extract_instruments(outcomes=e, p1=5e-8),silent = T)
    if (exists("exp_dat")==TRUE){
      exp_dat$samplesize.exposure=115078
      if (length(exp_dat$SNP)>0) {
        exp_dat$exposure=sapply(strsplit(exp_dat$exposure,fixed = TRUE, split= " ||"),"[",1)
        outcome_dat <-extract_outcome_data(
          snps = exp_dat$SNP,
          outcomes=n)
        if (length(outcome_dat$SNP)>0) {
          outcome_dat$outcome=sapply(strsplit(outcome_dat$outcome,fixed = TRUE, split= " ||"),"[",1)
          #My_MR(exp_dat,outcome_dat)
          rm(dat)
          try(dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat), silent=TRUE)
          dat=dat[is.na(dat$beta.outcome)==FALSE ,]
          
          ############# Change
          dat=dat[dat$mr_keep==TRUE,]
          
          if (length(which(dat$mr_keep=='TRUE'))>2) {
            #MendelianRandomization R Package
            MR_IVW_MRE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "random")
            MR_IVW_FE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "fixed")
            MR_Egger = MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                 by = dat$beta.outcome, byse = dat$se.outcome))
            MR_W_Med = MendelianRandomization::mr_median(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                  by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            MR_W_Mod = mr_mbe(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                       by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            
            dat_steiger = steiger_filtering(dat)
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.exposure")] <- "r.exposure"
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.outcome")] <- "r.outcome"
            
            dat_steiger_direct = directionality_test(dat_steiger)
            
            dat$len_SNP=length(dat$SNP)
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5082560/
            try(dat$R_Square<- 2*(dat$beta.exposure^2)*dat$eaf.exposure*(1-dat$eaf.exposure), silent=TRUE)
            try(dat$Total_R_Square <- sum(dat$R_Square), silent=TRUE)
            #try(dat$F_stat_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            try(dat$Total_R_Square <- dat$Total_R_Square/(dat$SD^2), silent=TRUE)
            try(dat$F_stat <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            #try(dat$F_stat_sim_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_stat_sim <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_statistic_pre <- (dat$beta.exposure/dat$se.exposure)^2, silent=TRUE)
            try(dat$F_statistic <- sum(dat$F_statistic_pre)/(dat$len_SNP), silent=TRUE)
            try(F_stat <- data.frame(dat$id.exposure[1],dat$len_SNP[1],dat$samplesize.exposure[1],dat$Total_R_Square[1],dat$F_stat[1],dat$F_stat_sim[1],dat$F_statistic[1]), silent=TRUE)
            
            if (length(which(dat$mr_keep=='TRUE'))>3) {
              PressoObject=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                                     SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, 
                                     NbDistribution = max(1000,length(dat$SNP)/0.05),  SignifThreshold = 0.05)
              PRESSO = PressoObject$`Main MR results`
              PRESSO_raw = PRESSO[PRESSO$"MR Analysis"=='Raw',]
              PRESSO_corrected = PRESSO[PRESSO$"MR Analysis"=='Outlier-corrected',]
              PRESSO_GLOBAL = PressoObject$`MR-PRESSO results`$`Global Test`
              pval_PRESSO_Global = PRESSO_GLOBAL$Pvalue
              names(PRESSO_raw)[names(PRESSO_raw) == "Causal Estimate"] <- "b_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "Sd"] <- "se_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "P-value"] <- "pval_PRESSO_raw"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Causal Estimate"] <- "b_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Sd"] <- "se_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "P-value"] <- "pval_PRESSO_corrected"
              
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   PRESSO_raw$b_PRESSO_raw, PRESSO_raw$se_PRESSO_raw, PRESSO_raw$pval_PRESSO_raw, 
                                   PRESSO_corrected$b_PRESSO_corrected, PRESSO_corrected$se_PRESSO_corrected, 
                                   PRESSO_corrected$pval_PRESSO_corrected, pval_PRESSO_Global)
              
              write.table(final_res, file= paste(Pathway_out,"MenRan_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
              
            } else {
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   NA,NA,NA,NA,NA,NA,NA)
              
              write.table(final_res, file= paste(Pathway_out,"MenRan_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
            }
          }
        }
      }
    }
  }
}
    
# ml R/4.2.1-foss-2020b
#################################################1e-6
library(TwoSampleMR)
library(dplyr)
'%ni%' <- Negate('%in%')
library(MendelianRandomization)
library(MRPRESSO)

ao <- available_outcomes()

Pathway_out="/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/MR_result/result_met_MS_06_22_22/"

####### met-a
data_met_a=ao
data_met_a$new_name_1=sapply(strsplit(data_met_a$id, split= "-"),"[",1)
data_met_a$new_name_2=sapply(strsplit(data_met_a$id, split= "-"),"[",2)
data_met_a=data_met_a[data_met_a$new_name_1=="met",]
data_met_a=data_met_a[data_met_a$new_name_2=="a",]

####### met-c
data_met_c=ao
data_met_c$new_name_1=sapply(strsplit(data_met_c$id, split= "-"),"[",1)
data_met_c$new_name_2=sapply(strsplit(data_met_c$id, split= "-"),"[",2)
data_met_c=data_met_c[data_met_c$new_name_1=="met",]
data_met_c=data_met_c[data_met_c$new_name_2=="c",]

####### met-d
data_met_d=ao
data_met_d$new_name_1=sapply(strsplit(data_met_d$id, split= "-"),"[",1)
data_met_d$new_name_2=sapply(strsplit(data_met_d$id, split= "-"),"[",2)
data_met_d=data_met_d[data_met_d$new_name_1=="met",]
data_met_d=data_met_d[data_met_d$new_name_2=="d",]

#b_Wald, se_Wald, pval_Wald,
final_res=data.frame("id.exposure","exposure","id.outcome","outcome",
                     "b_IVW_MRE", "se_IVW_MRE", "pval_IVW_MRE", "b_IVW_FE", "se_IVW_FE", "pval_IVW_FE",
                     "b_Egger", "se_Egger", "pval_Egger",
                     "Egger_intercept", "pval_intercept", "Het_IVW_pval", "Het_Egger_pval",
                     "b_W_Med", "se_W_Med", "pval_W_Med", "b_W_Mod", "se_W_Mod", "pval_W_Mod",
                     "nsnps","Total_R_Square", "F_stat", "F_stat_sim", "F_statistic",
                     "snp_r2.exposure","snp_r2.outcome","causal_direction","steiger_pval",
                     "b_PRESSO_raw", "se_PRESSO_raw", "pval_PRESSO_raw", 
                     "b_PRESSO_corrected", "se_PRESSO_corrected", "pval_PRESSO_corrected", "pval_PRESSO_Global")

# write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
#             row.names = F, quote = FALSE, na = "-",sep='\t')
write.table(final_res, file= paste(Pathway_out,"MenRan_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')

all_data_met <- c(data_met_a$id, data_met_c$id)
# all_data_met <- all_data_met[-21]
# all_data_met <- all_data_met[-20]
# all_data_met <- all_data_met[-19]
# 
# print(all_data_met)
options(scipen = 5)

for (e in all_data_met) {
  for (n in c("ieu-b-18")) {
    # > min(exp_dat$pval.exposure)
    # [1] 0.000001373
    # > 0.000001373*1e6
    # [1] 1.373
    if (e != "met-a-372") {
      try(exp_dat <- extract_instruments(outcomes=e, p1=1e-6),silent = T)
      if (exists("exp_dat")==TRUE){
        if (length(exp_dat$SNP)>0) {
          exp_dat$exposure=sapply(strsplit(exp_dat$exposure,fixed = TRUE, split= " ||"),"[",1)
          outcome_dat <-extract_outcome_data(
            snps = exp_dat$SNP,
            outcomes=n)
          if (length(outcome_dat$SNP)>0) {
            outcome_dat$outcome=sapply(strsplit(outcome_dat$outcome,fixed = TRUE, split= " ||"),"[",1)
            #My_MR(exp_dat,outcome_dat)
            rm(dat)
            try(dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat), silent=TRUE)
            dat=dat[is.na(dat$beta.outcome)==FALSE ,]
            
            ############# Change
            dat=dat[dat$mr_keep==TRUE,]
            
            if (length(which(dat$mr_keep=='TRUE'))>2) {
              #MendelianRandomization R Package
              MR_IVW_MRE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "random")
              MR_IVW_FE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "fixed")
              MR_Egger = MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                   by = dat$beta.outcome, byse = dat$se.outcome))
              MR_W_Med = MendelianRandomization::mr_median(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                    by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
              MR_W_Mod = mr_mbe(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                         by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
              
              dat_steiger = steiger_filtering(dat)
              colnames(dat_steiger)[which(names(dat_steiger) == "rsq.exposure")] <- "r.exposure"
              colnames(dat_steiger)[which(names(dat_steiger) == "rsq.outcome")] <- "r.outcome"
              
              dat_steiger_direct = directionality_test(dat_steiger)
              
              dat$len_SNP=length(dat$SNP)
              #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5082560/
              try(dat$R_Square<- 2*(dat$beta.exposure^2)*dat$eaf.exposure*(1-dat$eaf.exposure), silent=TRUE)
              try(dat$Total_R_Square <- sum(dat$R_Square), silent=TRUE)
              #try(dat$F_stat_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
              try(dat$Total_R_Square <- dat$Total_R_Square/(dat$SD^2), silent=TRUE)
              try(dat$F_stat <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
              #try(dat$F_stat_sim_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
              try(dat$F_stat_sim <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
              try(dat$F_statistic_pre <- (dat$beta.exposure/dat$se.exposure)^2, silent=TRUE)
              try(dat$F_statistic <- sum(dat$F_statistic_pre)/(dat$len_SNP), silent=TRUE)
              try(F_stat <- data.frame(dat$id.exposure[1],dat$len_SNP[1],dat$samplesize.exposure[1],dat$Total_R_Square[1],dat$F_stat[1],dat$F_stat_sim[1],dat$F_statistic[1]), silent=TRUE)
              
              if (length(which(dat$mr_keep=='TRUE'))>3) {
                PressoObject=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                                       SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, 
                                       NbDistribution = max(1000,length(dat$SNP)/0.05),  SignifThreshold = 0.05)
                PRESSO = PressoObject$`Main MR results`
                PRESSO_raw = PRESSO[PRESSO$"MR Analysis"=='Raw',]
                PRESSO_corrected = PRESSO[PRESSO$"MR Analysis"=='Outlier-corrected',]
                PRESSO_GLOBAL = PressoObject$`MR-PRESSO results`$`Global Test`
                pval_PRESSO_Global = PRESSO_GLOBAL$Pvalue
                names(PRESSO_raw)[names(PRESSO_raw) == "Causal Estimate"] <- "b_PRESSO_raw"
                names(PRESSO_raw)[names(PRESSO_raw) == "Sd"] <- "se_PRESSO_raw"
                names(PRESSO_raw)[names(PRESSO_raw) == "P-value"] <- "pval_PRESSO_raw"
                names(PRESSO_corrected)[names(PRESSO_corrected) == "Causal Estimate"] <- "b_PRESSO_corrected"
                names(PRESSO_corrected)[names(PRESSO_corrected) == "Sd"] <- "se_PRESSO_corrected"
                names(PRESSO_corrected)[names(PRESSO_corrected) == "P-value"] <- "pval_PRESSO_corrected"
                
                final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                     MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                     MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                     MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                     MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                     MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                     MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                     MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                     MR_IVW_FE@SNPs,
                                     dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                     dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                     dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                     PRESSO_raw$b_PRESSO_raw, PRESSO_raw$se_PRESSO_raw, PRESSO_raw$pval_PRESSO_raw, 
                                     PRESSO_corrected$b_PRESSO_corrected, PRESSO_corrected$se_PRESSO_corrected, 
                                     PRESSO_corrected$pval_PRESSO_corrected, pval_PRESSO_Global)
                
                # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                #             row.names = F, quote = FALSE, na = "-",sep='\t')
                write.table(final_res, file= paste(Pathway_out,"MenRan_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
                            row.names = F, quote = FALSE, na = "-",sep='\t')
                
              } else {
                final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                     MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                     MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                     MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                     MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                     MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                     MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                     MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                     MR_IVW_FE@SNPs,
                                     dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                     dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                     dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                     NA,NA,NA,NA,NA,NA,NA)
                
                # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
                #             row.names = F, quote = FALSE, na = "-",sep='\t')
                write.table(final_res, file= paste(Pathway_out,"MenRan_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
                            row.names = F, quote = FALSE, na = "-",sep='\t')
              }
            }
          }
        }
      }
    }
  }
}

all_data_met <- c(data_met_d$id)
# all_data_met <- all_data_met[-21]
# all_data_met <- all_data_met[-20]
# all_data_met <- all_data_met[-19]
# 
# print(all_data_met)

for (e in all_data_met) {
  for (n in c("ieu-b-18")) {
    try(exp_dat <- extract_instruments(outcomes=e, p1=1e-6),silent = T)
    if (exists("exp_dat")==TRUE){
      exp_dat$samplesize.exposure=115078
      if (length(exp_dat$SNP)>0) {
        exp_dat$exposure=sapply(strsplit(exp_dat$exposure,fixed = TRUE, split= " ||"),"[",1)
        outcome_dat <-extract_outcome_data(
          snps = exp_dat$SNP,
          outcomes=n)
        if (length(outcome_dat$SNP)>0) {
          outcome_dat$outcome=sapply(strsplit(outcome_dat$outcome,fixed = TRUE, split= " ||"),"[",1)
          #My_MR(exp_dat,outcome_dat)
          rm(dat)
          try(dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat), silent=TRUE)
          dat=dat[is.na(dat$beta.outcome)==FALSE ,]
          
          ############# Change
          dat=dat[dat$mr_keep==TRUE,]
          
          if (length(which(dat$mr_keep=='TRUE'))>2) {
            #MendelianRandomization R Package
            MR_IVW_MRE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "random")
            MR_IVW_FE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "fixed")
            MR_Egger = MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                 by = dat$beta.outcome, byse = dat$se.outcome))
            MR_W_Med = MendelianRandomization::mr_median(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                  by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            MR_W_Mod = mr_mbe(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                       by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
            
            dat_steiger = steiger_filtering(dat)
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.exposure")] <- "r.exposure"
            colnames(dat_steiger)[which(names(dat_steiger) == "rsq.outcome")] <- "r.outcome"
            
            dat_steiger_direct = directionality_test(dat_steiger)
            
            dat$len_SNP=length(dat$SNP)
            #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5082560/
            try(dat$R_Square<- 2*(dat$beta.exposure^2)*dat$eaf.exposure*(1-dat$eaf.exposure), silent=TRUE)
            try(dat$Total_R_Square <- sum(dat$R_Square), silent=TRUE)
            #try(dat$F_stat_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            try(dat$Total_R_Square <- dat$Total_R_Square/(dat$SD^2), silent=TRUE)
            try(dat$F_stat <- dat$Total_R_Square*(dat$samplesize.exposure-dat$len_SNP-1)/((1-dat$Total_R_Square)*(dat$len_SNP)), silent=TRUE)
            #try(dat$F_stat_sim_no_SD <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_stat_sim <- dat$Total_R_Square*(dat$samplesize.exposure)/(dat$len_SNP), silent=TRUE)
            try(dat$F_statistic_pre <- (dat$beta.exposure/dat$se.exposure)^2, silent=TRUE)
            try(dat$F_statistic <- sum(dat$F_statistic_pre)/(dat$len_SNP), silent=TRUE)
            try(F_stat <- data.frame(dat$id.exposure[1],dat$len_SNP[1],dat$samplesize.exposure[1],dat$Total_R_Square[1],dat$F_stat[1],dat$F_stat_sim[1],dat$F_statistic[1]), silent=TRUE)
            
            if (length(which(dat$mr_keep=='TRUE'))>3) {
              PressoObject=mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", 
                                     SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, 
                                     NbDistribution = max(1000,length(dat$SNP)/0.05),  SignifThreshold = 0.05)
              PRESSO = PressoObject$`Main MR results`
              PRESSO_raw = PRESSO[PRESSO$"MR Analysis"=='Raw',]
              PRESSO_corrected = PRESSO[PRESSO$"MR Analysis"=='Outlier-corrected',]
              PRESSO_GLOBAL = PressoObject$`MR-PRESSO results`$`Global Test`
              pval_PRESSO_Global = PRESSO_GLOBAL$Pvalue
              names(PRESSO_raw)[names(PRESSO_raw) == "Causal Estimate"] <- "b_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "Sd"] <- "se_PRESSO_raw"
              names(PRESSO_raw)[names(PRESSO_raw) == "P-value"] <- "pval_PRESSO_raw"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Causal Estimate"] <- "b_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "Sd"] <- "se_PRESSO_corrected"
              names(PRESSO_corrected)[names(PRESSO_corrected) == "P-value"] <- "pval_PRESSO_corrected"
              
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   PRESSO_raw$b_PRESSO_raw, PRESSO_raw$se_PRESSO_raw, PRESSO_raw$pval_PRESSO_raw, 
                                   PRESSO_corrected$b_PRESSO_corrected, PRESSO_corrected$se_PRESSO_corrected, 
                                   PRESSO_corrected$pval_PRESSO_corrected, pval_PRESSO_Global)
              
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
              write.table(final_res, file= paste(Pathway_out,"MenRan_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
              
            } else {
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs,
                                   dat$Total_R_Square[1], dat$F_stat[1], dat$F_stat_sim[1], dat$F_statistic[1],
                                   dat_steiger_direct$snp_r2.exposure,dat_steiger_direct$snp_r2.outcome,
                                   dat_steiger_direct$correct_causal_direction,dat_steiger_direct$steiger_pval,
                                   NA,NA,NA,NA,NA,NA,NA)
              
              # write.table(final_res, file= paste(Pathway_out,"MendelianRandomization_5e8.txt", sep=""), col.names = FALSE, append = TRUE,
              #             row.names = F, quote = FALSE, na = "-",sep='\t')
              write.table(final_res, file= paste(Pathway_out,"MenRan_1e6.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
            }
          }
        }
      }
    }
  }
}
  
############## Name

write.table(data_met_a, file= paste(Pathway_out,"Name.txt", sep=""), col.names = T, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')
write.table(data_met_c, file= paste(Pathway_out,"Name.txt", sep=""), col.names = F, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')
write.table(data_met_d, file= paste(Pathway_out,"Name.txt", sep=""), col.names = F, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')


############## Combine
name <- read.table(paste(Pathway_out,"Name.txt", sep=""), header=T, sep="\t")
name$id.exposure=name$id
name_id<-name%>%select(id.exposure, trait)

MenRan_1e6 <- read.table(paste(Pathway_out,"MenRan_1e6.txt", sep=""), header=T, sep="\t")

MenRan_1e6<- MenRan_1e6 %>% left_join(name_id, by= "id.exposure")
MenRan_1e6$exposure=MenRan_1e6$trait
MenRan_1e6 <- subset(MenRan_1e6, select = -trait)
write.table(MenRan_1e6, file= paste(Pathway_out,"MenRan_1e6_with_names.txt", sep=""), col.names = T, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')

MenRan_5e8 <- read.table(paste(Pathway_out,"MenRan_5e8.txt", sep=""), header=T, sep="\t")
MenRan_1e6 <- read.table(paste(Pathway_out,"MenRan_1e6_with_names.txt", sep=""), header=T, sep="\t")

MenRan_all<- MenRan_5e8 %>% left_join(MenRan_1e6, by= "id.exposure")

MenRan_all[MenRan_all$pval_IVW_MRE.x<0.05,]
MenRan_all[MenRan_all$pval_IVW_MRE.y<0.05,]

write.table(MenRan_all, file= paste(Pathway_out,"MenRan_all.txt", sep=""), col.names = T, append = TRUE,
            row.names = F, quote = FALSE, na = "-",sep='\t')


############## Compare
Yitang_5e8 <- read.table(paste("/Users/yelab/Cooperator/Metabolites_multiple_sclerosis_MR/","MenRan_5e8.txt", sep=""), header=T, sep="\t")

Angela_5e8 <- read.table(paste("/Users/yelab/Cooperator/Metabolites_multiple_sclerosis_MR/","TwoSampleMR_vs_MendelianRandomization_2.txt", sep=""), header=T, sep="\t")

Yitang_5e8$id.exposure=as.character(Yitang_5e8$id.exposure)

Yitang_5e8[duplicated(Yitang_5e8$id.exposure),]$id.exposure

data_met_a[duplicated(data_met_a$id),]$id
data_met_c[duplicated(data_met_c$id),]$id
data_met_d[duplicated(data_met_d$id),]$id

data_met=rbind(data_met_a,data_met_c,data_met_d)
data_met[duplicated(data_met$id),]$id

Yitang_1e6 <- read.table(paste("/Users/yelab/Cooperator/Metabolites_multiple_sclerosis_MR/","MenRan_1e6_with_names.txt", sep=""), header=T, sep="\t")

Angela_1e6 <- read.table(paste("/Users/yelab/Cooperator/Metabolites_multiple_sclerosis_MR/","TwoSampleMR_vs_MendelianRandomization_pval_change_2.txt", sep=""), header=T, sep="\t")

Yitang_1e6$id.exposure=as.character(Yitang_1e6$id.exposure)

Yitang_1e6[duplicated(Yitang_1e6$id.exposure),]$id.exposure

Angela_1e6$id.exposure=as.character(Angela_1e6$id.exposure)

Angela_1e6[duplicated(Angela_1e6$id.exposure),]$id.exposure


Comparison_5e8<- Yitang_5e8 %>% full_join(Angela_5e8, by= "id.exposure")

Comparison_5e8$compare_pval_IVW_MRE= Comparison_5e8$pval_IVW_MRE.x-Comparison_5e8$pval_IVW_MRE.y
summary(Comparison_5e8$compare_pval_IVW_MRE)

Comparison_5e8$pval_Egger= Comparison_5e8$pval_Egger.x-Comparison_5e8$pval_Egger.y
summary(Comparison_5e8$pval_Egger)


