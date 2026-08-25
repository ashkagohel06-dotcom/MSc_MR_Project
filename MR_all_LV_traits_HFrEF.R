library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(openxlsx)

setwd("C:/Users/HP/Documents/msc_project")

# ==============================================================================
# DEFINE ALL LV TRAITS
# ==============================================================================
lv_traits <- c(
  "LVEDV",
  "LVEDV_BSA",
  "LVEF",
  "LVESV",
  "LVESV_BSA",
  "LVGFI",
  "LVM",
  "LVM_BSA",
  "LVMCF",
  "LVMVR",
  "LVSV",
  "LVSV_BSA"
)

# Create empty results dataframe
all_results <- list()

# ==============================================================================
# LOOP THROUGH EACH TRAIT
# ==============================================================================
for (trait in lv_traits) {
  
  message(paste("Processing trait:", trait))
  
  # ==============================================================================
  # LOAD GWAS DATA
  # ==============================================================================
  gwas_path <- paste0(
    "C:/Users/HP/Documents/msc_project/input/Data/Data/LV_GWAS_Data/",
    trait,
    "_38_37_rsids.txt"
  )
  
  gwas <- fread(gwas_path)
  
  # ==============================================================================
  # FILTER SIGNIFICANT SNPs
  # ==============================================================================
  sig_snps <- gwas %>% 
    filter(P < 5e-8) %>%
    select(rsid_1kg, P)
  
  fwrite(sig_snps, 'temp_to_clump.txt', sep='\t', col.names = TRUE)
  
  # ==============================================================================
  # PLINK CLUMPING
  # ==============================================================================
  plink_cmd <- paste(
    "C:/Users/HP/Downloads/plink_win64_20250819/plink.exe",
    "--bfile C:/Users/HP/Documents/msc_project/input/Data/Data/1000Genomes_Project_Data/EUR_phase3_autosomes",
    "--clump temp_to_clump.txt",
    "--clump-p1 5e-8",
    "--clump-p2 5e-8",
    "--clump-r2 0.001",
    "--clump-kb 10000",
    "--clump-field P",
    "--clump-snp-field rsid_1kg",
    "--out clump_result"
  )
  
  system(plink_cmd)
  cat(readLines("clump_result.log"), sep = "\n")
  
  # ==============================================================================
  # EXTRACT INDEPENDENT SNPs
  # ==============================================================================
  clumped_list <- fread("clump_result.clumped")
  final_instruments <- gwas %>% 
    filter(rsid_1kg %in% clumped_list$SNP)
  
  fwrite(final_instruments, "clumped_lead_snps.csv")
  
  # ==============================================================================
  # F-STATISTICS
  # ==============================================================================
  f_stats_table <- data.table(
    rsid = final_instruments$rsid_1kg,
    f_stat = (final_instruments$BETA^2) / (final_instruments$SE^2)
  )
  
  f_mean   <- mean(f_stats_table$f_stat, na.rm = TRUE)
  f_median <- median(f_stats_table$f_stat, na.rm = TRUE)
  
  message("Individual F-statistics:")
  print(f_stats_table)
  
  message("\nSummary Statistics:")
  cat(sprintf("Mean F-statistic:   %.2f\n", f_mean))
  cat(sprintf("Median F-statistic: %.2f\n", f_median))
  
  if (f_mean < 10) {
    warning("CAUTION: Instruments may be weak")
  } else {
    message("PASS: Instruments are strong")
  }
  
  # ==============================================================================
  # MENDELIAN RANDOMISATION
  # ==============================================================================
  
  exp_raw <- fread("clumped_lead_snps.csv") %>%
    as.data.frame()
  
  out_full <- fread("C:/Users/HP/Documents/msc_project/input/HF_Data/HF_Data/HFrEF_GCST90726620.tsv") %>%
    as.data.frame()
  
  # FORMAT DATA
  exp_dat <- format_data(
    exp_raw,
    type = "exposure",
    snp_col = "rsid_1kg",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "MAF",
    pval_col = "P",
    phenotype_col = trait
  )
  
  out_dat <- format_data(
    out_full,
    type = "outcome",
    snps = exp_dat$SNP,
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    phenotype_col = "HrEFF"
  )
  
  # HARMONISE
  dat <- harmonise_data(exp_dat, out_dat, action = 2)
  n_snps <- nrow(dat)
  
  # MR INPUT
  mr_in <- mr_input(
    bx = dat$beta.exposure,
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome,
    snps = dat$SNP
  )
  
  # IVW + HETEROGENEITY
  ivw <- tryCatch(
    mr_ivw(mr_in, robust = TRUE, penalized = TRUE),
    error = function(e) NULL
  )
  
  het <- tryCatch(
    TwoSampleMR::mr_heterogeneity(dat),
    error = function(e) NULL
  )
  
  # ==============================================================================
  # FINAL RESULTS TABLE
  # ==============================================================================
  results <- tibble(
    Phenotype = trait,
    SNP_count = n_snps,
    IVW_beta = if(!is.null(ivw)) ivw@Estimate else NA,
    IVW_CI_low = if(!is.null(ivw)) ivw@CILower else NA,
    IVW_CI_high = if(!is.null(ivw)) ivw@CIUpper else NA,
    IVW_p = if(!is.null(ivw)) ivw@Pvalue else NA,
    Q_p = if(!is.null(het) && nrow(het) > 0) het$Q_pval[1] else NA
  )
  
  # CALCULATE OR
  results <- results %>%
    mutate(
      IVW_OR = exp(IVW_beta),
      IVW_OR_LCI = exp(IVW_CI_low),
      IVW_OR_UCI = exp(IVW_CI_high)
    ) %>%
    select(
      Phenotype,
      SNP_count,
      IVW_beta,
      IVW_p,
      IVW_OR,
      IVW_OR_LCI,
      IVW_OR_UCI,
      Q_p
    )
  
  print(results)
  
  # STORE RESULTS
  all_results[[trait]] <- results
}

# ==============================================================================
# COMBINE ALL RESULTS
# ==============================================================================
final_output <- bind_rows(all_results)

print(final_output)

# ==============================================================================
# SAVE FINAL EXCEL FILE
# ==============================================================================
write.xlsx(
  final_output,
  "C:/Users/HP/Documents/msc_project/output/all_LV_traits_IVW_results_HFrEF.xlsx",
  overwrite = TRUE
)

message("All traits processed successfully!")