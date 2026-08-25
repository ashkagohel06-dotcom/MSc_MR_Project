library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(openxlsx)

setwd("C:/Users/HP/Documents/msc_project")

# ==============================================================================
# DEFINE ALL RV TRAITS
# ==============================================================================
rv_traits <- c(
  "RV_LV_ratio",
  "RVEDV",
  "RVEDV_BSA",
  "RVEF",
  "RVESV",
  "RVESV_BSA",
  "RVSV",
  "RVSV_BSA"
)

# HF outcomes (TSV files, same as LV work)
hf_outcomes <- list(
  allcause_HF = "input/HF_Data/HF_Data/allcause_HF_GCST90726617.tsv",
  HFpEF      = "input/HF_Data/HF_Data/HFpEF_GCST90726619.tsv",
  HFrEF      = "input/HF_Data/HF_Data/HFrEF_GCST90726620.tsv"
)

all_results <- list()

# ==============================================================================
# LOOP THROUGH EACH RV TRAIT × HF OUTCOME
# ==============================================================================
for (trait in rv_traits) {
  for (hf_name in names(hf_outcomes)) {
    
    message(paste("Processing:", trait, "→", hf_name))
    
    # ==============================================================================
    # LOAD RV GWAS DATA
    # ==============================================================================
    gwas_path <- paste0(
      "C:/Users/HP/Downloads/RV_GWAS_Data/",
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
    
    fwrite(sig_snps, "temp_to_clump.txt", sep = "\t", col.names = TRUE)
    
    # ==============================================================================
    # PLINK CLUMPING (same as LV)
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
    
    out_full <- fread(hf_outcomes[[hf_name]]) %>%
      as.data.frame()
    
    # FORMAT EXPOSURE (RV trait)
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
    
    # FORMAT OUTCOME (HF subtype)
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
      phenotype_col = hf_name
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
      RV_Phenotype = trait,
      Outcome      = hf_name,
      SNP_count    = n_snps,
      IVW_beta     = if (!is.null(ivw)) ivw@Estimate else NA,
      IVW_CI_low   = if (!is.null(ivw)) ivw@CILower else NA,
      IVW_CI_high  = if (!is.null(ivw)) ivw@CIUpper else NA,
      IVW_p        = if (!is.null(ivw)) ivw@Pvalue else NA,
      Q_p          = if (!is.null(het) && nrow(het) > 0) het$Q_pval[1] else NA
    ) %>%
      mutate(
        IVW_OR      = exp(IVW_beta),
        IVW_OR_LCI  = exp(IVW_CI_low),
        IVW_OR_UCI  = exp(IVW_CI_high)
      ) %>%
      select(
        RV_Phenotype,
        Outcome,
        SNP_count,
        IVW_beta,
        IVW_p,
        IVW_OR,
        IVW_OR_LCI,
        IVW_OR_UCI,
        Q_p
      )
    
    print(results)
    
    all_results[[paste(trait, hf_name, sep = "_")]] <- results
  }
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
  "C:/Users/HP/Documents/msc_project/output/all_RV_traits_IVW_results_allHF.xlsx",
  overwrite = TRUE
)

message("All RV traits × HF outcomes processed successfully!")
