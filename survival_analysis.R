library(maftools)
library(tidyverse)
library(stringr)
setwd("E:/share_hdd/research/T4GA/MAF/")


raw_files <- list.files("E:/share_hdd/research/T4GA/MAF/maf", "maf.gz", full.names = T, recursive = T)
clinical_df <- read.delim("E:/share_hdd/research/T4GA/MAF/clinical/clinical.tsv")


clin <- clinical_df %>% 
  dplyr::select(c(
    "case_submitter_id", "days_to_last_follow_up", "vital_status",
    "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
    "ajcc_pathologic_stage", "gender"
  )) %>%
  `names<-`(c("Tumor_Sample_Barcode", "time", "status", "T", "N", "M", "stage", "gender")) %>% distinct()

#merged_maf <- merge_mafs(raw_files,verbose = T)


merged_maf <- read.maf("merged.maf", removeDuplicatedVariants=F)


tmp_df  <- merged_maf@maf.silent[merged_maf@maf.silent[["Hugo_Symbol"]]=="BRCA2",]
tmp_df$Tumor_Sample_Barcode_bk <- tmp_df$Tumor_Sample_Barcode
tmp_df$Tumor_Sample_Barcode <- str_sub(tmp_df$Tumor_Sample_Barcode, start = 1, end = 12)

tmp_df <- left_join(tmp_df,clin,by="Tumor_Sample_Barcode") %>% 
  subset(stage %in% c("Stage IIA","Stage IIB") & time != "'--") %>% 
  mutate(label=if_else(dbSNP_RS=="rs15869", 1, 0),
         time = as.numeric(.$time))

clin_df <- tmp_df %>% dplyr::select(c("Tumor_Sample_Barcode", "time", "status", "label", "Tumor_Sample_Barcode_bk")) %>% distinct()


library(survival)
library(survminer)

data.surv <- clin_df %>% filter(!is.na(time)) %>%
  mutate(
    time = time / 30, 
    status = if_else(status == "Alive", 0, 1)
  )

km_fit <- survfit(Surv(time, status) ~ label, data = data.surv)

ggsurvplot(
  km_fit, data = data.surv,
  pval = TRUE, surv.median.line = "hv",
  legend.labs=c("w/o rs15869","w/ rs15869"),
  legend.title="Group",
  title="Overall survival",
  xlab = "Time(month)",
  risk.table = TRUE
)


# TMB analysis
maf_sub <- subsetMaf(merged_maf, tsb = clin_df$Tumor_Sample_Barcode_bk)

BRCA2_tmb <- tmb(maf_sub)

clin_df <- left_join(clin_df, BRCA2_tmb, by = c("Tumor_Sample_Barcode_bk" = "Tumor_Sample_Barcode"))

clin_df <- clin_df %>% arrange(total_perMB_log) %>% 
  mutate(SNP_Status = if_else(label==1,"w/ rs15869","w/o rs15869"),
         seq = seq(nrow(.)))

ggplot(clin_df,aes(x=seq, y=total_perMB_log, color=SNP_Status)) + 
  geom_point()+
  theme_classic()+
  geom_abline(slope=0,intercept = median(clin_df$total_perMB_log),color="red")+
  labs(title = "Mutation Burden") + 
  ylab("TMB/MB(log10)") + xlab("sample Median: 19.94/MB")


ggplot(clin_df,aes(x=SNP_Status, y=total_perMB_log, fill=SNP_Status)) + 
  geom_boxplot(width=0.4)+
  geom_point() +
  theme_classic()+
  labs(title = "Mutation Burden Boxplot") + 
  ylab("TMB/MB(log10)") + xlab("SNP status")

wilcox.test(x= clin_df$total_perMB[clin_df$label==1],y=clin_df$total_perMB[!clin_df$label==1])
