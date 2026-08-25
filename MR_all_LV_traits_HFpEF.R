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
  "LVMCF",# ==============================================================================
  # AUTOMATED TWO-SAMPLE MENDELIAN RANDOMISATION
  # LV TRAITS -> HF OUTCOME
  #
  # Analyses:
  #   1. IVW
  #   2. MR-Egger
  #   3. Weighted Median
  #   4. Steiger directionality
  #   5. F-statistics
  #   6. Cochran's Q heterogeneity
  #   7. MR-Egger intercept / pleiotropy
  #
  # This script is designed to be run separately for:
  #   - HFpEF
  #   - HFrEF
  #   - All-cause HF
  #
  # CHANGE ONLY THE OUTCOME SETTINGS BELOW WHEN COPYING THIS SCRIPT
  # ==============================================================================
  
  
  # ==============================================================================
  # 1. LOAD PACKAGES
  # ==============================================================================
  
  library(data.table)
  library(tidyverse)
  library(TwoSampleMR)
  library(openxlsx)
  
  
  # ==============================================================================
  # 2. WORKING DIRECTORY
  # ==============================================================================
  
  setwd("C:/Users/HP/Documents/msc_project")
  
  
  # ==============================================================================
  # 3. DEFINE LV TRAITS
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
  
  
  # ==============================================================================
  # 4. OUTCOME SETTINGS
  #
  # CHANGE THESE SETTINGS FOR EACH OF THE THREE SCRIPTS
  # ==============================================================================
  
  outcome_name <- "HFpEF"
  
  outcome_file <- paste0(
    "C:/Users/HP/Documents/msc_project/input/HF_Data/HF_Data/",
    "HFpEF_GCST90726619.tsv"
  )
  
  output_file <- paste0(
    "C:/Users/HP/Documents/msc_project/output/",
    "all_LV_traits_IVW_results_",
    outcome_name,
    ".xlsx"
  )
  
  
  # ==============================================================================
  # 5. PLINK SETTINGS
  # ==============================================================================
  
  plink_exe <- "C:/Users/HP/Downloads/plink_win64_20250819/plink.exe"
  
  plink_bfile <- paste0(
    "C:/Users/HP/Documents/msc_project/input/Data/Data/",
    "1000Genomes_Project_Data/EUR_phase3_autosomes"
  )
  
  
  # ==============================================================================
  # 6. INPUT DIRECTORY FOR LV GWAS DATA
  # ==============================================================================
  
  lv_gwas_directory <- paste0(
    "C:/Users/HP/Documents/msc_project/input/Data/Data/",
    "LV_GWAS_Data/"
  )
  
  
  # ==============================================================================
  # 7. LOAD HF OUTCOME DATA ONCE
  # ==============================================================================
  
  message("Loading outcome GWAS: ", outcome_name)
  
  out_full <- fread(outcome_file) %>%
    as.data.frame()
  
  message(
    "Outcome GWAS loaded successfully: ",
    nrow(out_full),
    " rows"
  )
  
  
  # ==============================================================================
  # 8. CREATE EMPTY RESULT LISTS
  # ==============================================================================
  
  all_results <- list()
  
  all_f_statistics <- list()
  
  all_harmonised_data <- list()
  
  
  # ==============================================================================
  # 9. LOOP THROUGH EACH LV TRAIT
  # ==============================================================================
  
  for (trait in lv_traits) {
    
    message("")
    message("==============================================================")
    message("Processing trait: ", trait)
    message("Outcome: ", outcome_name)
    message("==============================================================")
    
    
    # ============================================================================
    # 9.1 LOAD LV GWAS
    # ============================================================================
    
    gwas_path <- paste0(
      lv_gwas_directory,
      trait,
      "_38_37_rsids.txt"
    )
    
    if (!file.exists(gwas_path)) {
      
      warning(
        "GWAS file does not exist for trait: ",
        trait
      )
      
      next
    }
    
    gwas <- fread(gwas_path)
    
    
    # ============================================================================
    # 9.2 CHECK REQUIRED LV COLUMNS
    # ============================================================================
    
    required_exposure_columns <- c(
      "rsid_1kg",
      "P",
      "BETA",
      "SE",
      "ALLELE1",
      "ALLELE0",
      "MAF"
    )
    
    missing_exposure_columns <- setdiff(
      required_exposure_columns,
      names(gwas)
    )
    
    if (length(missing_exposure_columns) > 0) {
      
      warning(
        "Missing exposure columns for ",
        trait,
        ": ",
        paste(missing_exposure_columns, collapse = ", ")
      )
      
      next
    }
    
    
    # ============================================================================
    # 9.3 SELECT GENOME-WIDE SIGNIFICANT SNPs
    # ============================================================================
    
    sig_snps <- gwas %>%
      filter(
        !is.na(P),
        P < 5e-8
      ) %>%
      select(
        rsid_1kg,
        P
      ) %>%
      distinct(rsid_1kg, .keep_all = TRUE)
    
    
    n_sig <- nrow(sig_snps)
    
    message(
      "Genome-wide significant SNPs: ",
      n_sig
    )
    
    
    # ============================================================================
    # STOP IF NO SIGNIFICANT SNPs
    # ============================================================================
    
    if (n_sig == 0) {
      
      warning(
        "No genome-wide significant SNPs for ",
        trait
      )
      
      next
    }
    
    
    # ============================================================================
    # 9.4 TRAIT-SPECIFIC PLINK FILES
    # ============================================================================
    
    temp_file <- paste0(
      trait,
      "_temp_to_clump.txt"
    )
    
    clump_prefix <- paste0(
      trait,
      "_clump_result"
    )
    
    clumped_file <- paste0(
      clump_prefix,
      ".clumped"
    )
    
    
    # ============================================================================
    # 9.5 WRITE SNPs FOR PLINK
    # ============================================================================
    
    fwrite(
      sig_snps,
      temp_file,
      sep = "\t",
      col.names = TRUE
    )
    
    
    # ============================================================================
    # 9.6 PLINK CLUMPING
    # ============================================================================
    
    plink_cmd <- paste(
      shQuote(plink_exe),
      
      "--bfile",
      shQuote(plink_bfile),
      
      "--clump",
      shQuote(temp_file),
      
      "--clump-p1 5e-8",
      
      "--clump-p2 5e-8",
      
      "--clump-r2 0.001",
      
      "--clump-kb 10000",
      
      "--clump-field P",
      
      "--clump-snp-field rsid_1kg",
      
      "--out",
      shQuote(clump_prefix)
    )
    
    
    message("Running PLINK clumping...")
    
    system(plink_cmd)
    
    
    # ============================================================================
    # 9.7 CHECK WHETHER PLINK PRODUCED OUTPUT
    # ============================================================================
    
    if (!file.exists(clumped_file)) {
      
      warning(
        "No PLINK clumped file produced for ",
        trait
      )
      
      next
    }
    
    
    # ============================================================================
    # 9.8 READ CLUMPED SNPs
    # ============================================================================
    
    clumped_list <- fread(
      clumped_file
    )
    
    
    n_clumped <- nrow(clumped_list)
    
    message(
      "Independent SNPs after clumping: ",
      n_clumped
    )
    
    
    if (n_clumped == 0) {
      
      warning(
        "No independent SNPs after clumping for ",
        trait
      )
      
      next
    }
    
    
    # ============================================================================
    # 9.9 EXTRACT INDEPENDENT SNPs FROM ORIGINAL GWAS
    # ============================================================================
    
    final_instruments <- gwas %>%
      filter(
        rsid_1kg %in% clumped_list$SNP
      ) %>%
      distinct(
        rsid_1kg,
        .keep_all = TRUE
      )
    
    
    # ============================================================================
    # 9.10 SAVE TRAIT-SPECIFIC INSTRUMENTS
    # ============================================================================
    
    fwrite(
      final_instruments,
      paste0(
        trait,
        "_clumped_lead_snps.csv"
      )
    )
    
    
    # ============================================================================
    # 9.11 F-STATISTICS
    # ============================================================================
    
    f_stats_table <- data.table(
      Phenotype = trait,
      
      rsid = final_instruments$rsid_1kg,
      
      BETA = final_instruments$BETA,
      
      SE = final_instruments$SE,
      
      P = final_instruments$P,
      
      F_statistic =
        (final_instruments$BETA^2) /
        (final_instruments$SE^2)
    )
    
    
    # Remove invalid F-statistics
    
    f_stats_table <- f_stats_table %>%
      filter(
        is.finite(F_statistic),
        !is.na(F_statistic)
      )
    
    
    f_mean <- mean(
      f_stats_table$F_statistic,
      na.rm = TRUE
    )
    
    f_median <- median(
      f_stats_table$F_statistic,
      na.rm = TRUE
    )
    
    f_min <- min(
      f_stats_table$F_statistic,
      na.rm = TRUE
    )
    
    f_max <- max(
      f_stats_table$F_statistic,
      na.rm = TRUE
    )
    
    
    message(
      "Mean F-statistic: ",
      round(f_mean, 2)
    )
    
    message(
      "Median F-statistic: ",
      round(f_median, 2)
    )
    
    message(
      "Minimum F-statistic: ",
      round(f_min, 2)
    )
    
    
    if (f_min < 10) {
      
      warning(
        "At least one instrument has F < 10 for ",
        trait
      )
      
    } else {
      
      message(
        "PASS: All instruments have F >= 10"
      )
    }
    
    
    # Store F statistics
    
    all_f_statistics[[trait]] <- f_stats_table
    
    
    # ============================================================================
    # 9.12 FORMAT EXPOSURE DATA
    # ============================================================================
    
    exp_raw <- final_instruments %>%
      as.data.frame()
    
    
    exp_dat <- tryCatch(
      
      format_data(
        
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
      ),
      
      error = function(e) {
        
        warning(
          "Exposure formatting failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (is.null(exp_dat)) {
      next
    }
    
    
    # ============================================================================
    # 9.13 FORMAT OUTCOME DATA
    # ============================================================================
    
    out_dat <- tryCatch(
      
      format_data(
        
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
        
        phenotype_col = outcome_name
      ),
      
      error = function(e) {
        
        warning(
          "Outcome formatting failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (is.null(out_dat)) {
      next
    }
    
    
    # ============================================================================
    # 9.14 HARMONISE
    # ============================================================================
    
    dat <- tryCatch(
      
      harmonise_data(
        exp_dat,
        out_dat,
        action = 2
      ),
      
      error = function(e) {
        
        warning(
          "Harmonisation failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (is.null(dat)) {
      next
    }
    
    
    # ============================================================================
    # 9.15 KEEP ONLY SNPs THAT CAN BE USED FOR MR
    # ============================================================================
    
    dat <- dat %>%
      filter(
        mr_keep == TRUE
      )
    
    
    n_harmonised <- nrow(dat)
    
    
    message(
      "Harmonised SNPs available for MR: ",
      n_harmonised
    )
    
    
    # ============================================================================
    # STOP IF INSUFFICIENT SNPs
    # ============================================================================
    
    if (n_harmonised == 0) {
      
      warning(
        "No usable SNPs after harmonisation for ",
        trait
      )
      
      next
    }
    
    
    # Save harmonised data
    
    all_harmonised_data[[trait]] <- dat
    
    
    # ============================================================================
    # 9.16 MR ANALYSIS
    #
    # IVW
    # MR-Egger
    # Weighted Median
    # ============================================================================
    
    mr_results <- tryCatch(
      
      mr(
        dat,
        
        method_list = c(
          "mr_ivw",
          "mr_egger_regression",
          "mr_weighted_median"
        )
      ),
      
      error = function(e) {
        
        warning(
          "MR analysis failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (is.null(mr_results)) {
      next
    }
    
    
    # ============================================================================
    # 9.17 EXTRACT IVW RESULT
    # ============================================================================
    
    ivw_result <- mr_results %>%
      filter(
        method == "Inverse variance weighted"
      )
    
    
    if (nrow(ivw_result) == 0) {
      
      ivw_result <- tibble(
        b = NA_real_,
        se = NA_real_,
        pval = NA_real_
      )
      
    } else {
      
      ivw_result <- ivw_result[1, ]
    }
    
    
    # ============================================================================
    # 9.18 EXTRACT MR-EGGER RESULT
    # ============================================================================
    
    egger_result <- mr_results %>%
      filter(
        method == "MR Egger"
      )
    
    
    if (nrow(egger_result) == 0) {
      
      egger_result <- tibble(
        b = NA_real_,
        se = NA_real_,
        pval = NA_real_
      )
      
    } else {
      
      egger_result <- egger_result[1, ]
    }
    
    
    # ============================================================================
    # 9.19 EXTRACT WEIGHTED MEDIAN RESULT
    # ============================================================================
    
    weighted_median_result <- mr_results %>%
      filter(
        method == "Weighted median"
      )
    
    
    if (nrow(weighted_median_result) == 0) {
      
      weighted_median_result <- tibble(
        b = NA_real_,
        se = NA_real_,
        pval = NA_real_
      )
      
    } else {
      
      weighted_median_result <-
        weighted_median_result[1, ]
    }
    
    
    # ============================================================================
    # 9.20 HETEROGENEITY
    # ============================================================================
    
    heterogeneity <- tryCatch(
      
      mr_heterogeneity(
        dat,
        
        method_list = c(
          "mr_ivw",
          "mr_egger_regression"
        )
      ),
      
      error = function(e) {
        
        warning(
          "Heterogeneity analysis failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    # ============================================================================
    # 9.21 EXTRACT IVW Q
    # ============================================================================
    
    if (!is.null(heterogeneity)) {
      
      ivw_het <- heterogeneity %>%
        filter(
          method == "Inverse variance weighted"
        )
      
    } else {
      
      ivw_het <- NULL
    }
    
    
    if (
      !is.null(ivw_het) &&
      nrow(ivw_het) > 0
    ) {
      
      Q <- ivw_het$Q[1]
      
      Q_df <- ivw_het$Q_df[1]
      
      Q_p <- ivw_het$Q_pval[1]
      
    } else {
      
      Q <- NA_real_
      
      Q_df <- NA_real_
      
      Q_p <- NA_real_
    }
    
    
    # ============================================================================
    # 9.22 MR-EGGER INTERCEPT / PLEIOTROPY TEST
    # ============================================================================
    
    pleiotropy <- tryCatch(
      
      mr_pleiotropy_test(dat),
      
      error = function(e) {
        
        warning(
          "MR-Egger intercept test failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (
      !is.null(pleiotropy) &&
      nrow(pleiotropy) > 0
    ) {
      
      egger_intercept <- pleiotropy$egger_intercept[1]
      
      egger_intercept_se <- pleiotropy$se[1]
      
      egger_intercept_p <- pleiotropy$pval[1]
      
    } else {
      
      egger_intercept <- NA_real_
      
      egger_intercept_se <- NA_real_
      
      egger_intercept_p <- NA_real_
    }
    
    
    # ============================================================================
    # 9.23 STEIGER DIRECTIONALITY TEST
    # ============================================================================
    
    steiger <- tryCatch(
      
      directionality_test(dat),
      
      error = function(e) {
        
        warning(
          "Steiger test failed for ",
          trait,
          ": ",
          e$message
        )
        
        return(NULL)
      }
    )
    
    
    if (
      !is.null(steiger) &&
      nrow(steiger) > 0
    ) {
      
      steiger_r_exposure <- steiger$r.exposure[1]
      
      steiger_r_outcome <- steiger$r.outcome[1]
      
      steiger_p <- steiger$steiger_pval[1]
      
      steiger_direction <- steiger$correct_causal_direction[1]
      
    } else {
      
      steiger_r_exposure <- NA_real_
      
      steiger_r_outcome <- NA_real_
      
      steiger_p <- NA_real_
      
      steiger_direction <- NA
    }
    
    
    # ============================================================================
    # 9.24 CALCULATE CONFIDENCE INTERVALS
    # ============================================================================
    
    ivw_beta <- ivw_result$b
    
    ivw_se <- ivw_result$se
    
    ivw_p <- ivw_result$pval
    
    
    ivw_ci_low <- ivw_beta -
      1.96 * ivw_se
    
    ivw_ci_high <- ivw_beta +
      1.96 * ivw_se
    
    
    ivw_or <- exp(ivw_beta)
    
    ivw_or_lci <- exp(ivw_ci_low)
    
    ivw_or_uci <- exp(ivw_ci_high)
    
    
    # MR-Egger
    
    egger_beta <- egger_result$b
    
    egger_se <- egger_result$se
    
    egger_p <- egger_result$pval
    
    
    egger_ci_low <- egger_beta -
      1.96 * egger_se
    
    egger_ci_high <- egger_beta +
      1.96 * egger_se
    
    
    egger_or <- exp(egger_beta)
    
    egger_or_lci <- exp(egger_ci_low)
    
    egger_or_uci <- exp(egger_ci_high)
    
    
    # Weighted Median
    
    wm_beta <- weighted_median_result$b
    
    wm_se <- weighted_median_result$se
    
    wm_p <- weighted_median_result$pval
    
    
    wm_ci_low <- wm_beta -
      1.96 * wm_se
    
    wm_ci_high <- wm_beta +
      1.96 * wm_se
    
    
    wm_or <- exp(wm_beta)
    
    wm_or_lci <- exp(wm_ci_low)
    
    wm_or_uci <- exp(wm_ci_high)
    
    
    # ============================================================================
    # 9.25 FINAL SUMMARY RESULT
    # ============================================================================
    
    results <- tibble(
      
      Phenotype = trait,
      
      Outcome = outcome_name,
      
      Significant_SNPs = n_sig,
      
      Clumped_SNPs = n_clumped,
      
      Harmonised_SNPs = n_harmonised,
      
      
      # --------------------------------------------------------------------------
      # F STATISTICS
      # --------------------------------------------------------------------------
      
      Mean_F = f_mean,
      
      Median_F = f_median,
      
      Minimum_F = f_min,
      
      Maximum_F = f_max,
      
      
      # --------------------------------------------------------------------------
      # IVW
      # --------------------------------------------------------------------------
      
      IVW_beta = ivw_beta,
      
      IVW_SE = ivw_se,
      
      IVW_CI_low = ivw_ci_low,
      
      IVW_CI_high = ivw_ci_high,
      
      IVW_p = ivw_p,
      
      IVW_OR = ivw_or,
      
      IVW_OR_LCI = ivw_or_lci,
      
      IVW_OR_UCI = ivw_or_uci,
      
      
      # --------------------------------------------------------------------------
      # MR-EGGER
      # --------------------------------------------------------------------------
      
      Egger_beta = egger_beta,
      
      Egger_SE = egger_se,
      
      Egger_CI_low = egger_ci_low,
      
      Egger_CI_high = egger_ci_high,
      
      Egger_p = egger_p,
      
      Egger_OR = egger_or,
      
      Egger_OR_LCI = egger_or_lci,
      
      Egger_OR_UCI = egger_or_uci,
      
      
      # --------------------------------------------------------------------------
      # WEIGHTED MEDIAN
      # --------------------------------------------------------------------------
      
      WeightedMedian_beta = wm_beta,
      
      WeightedMedian_SE = wm_se,
      
      WeightedMedian_CI_low = wm_ci_low,
      
      WeightedMedian_CI_high = wm_ci_high,
      
      WeightedMedian_p = wm_p,
      
      WeightedMedian_OR = wm_or,
      
      WeightedMedian_OR_LCI = wm_or_lci,
      
      WeightedMedian_OR_UCI = wm_or_uci,
      
      
      # --------------------------------------------------------------------------
      # HETEROGENEITY
      # --------------------------------------------------------------------------
      
      IVW_Q = Q,
      
      IVW_Q_df = Q_df,
      
      IVW_Q_p = Q_p,
      
      
      # --------------------------------------------------------------------------
      # MR-EGGER PLEIOTROPY
      # --------------------------------------------------------------------------
      
      Egger_intercept = egger_intercept,
      
      Egger_intercept_SE = egger_intercept_se,
      
      Egger_intercept_p = egger_intercept_p,
      
      
      # --------------------------------------------------------------------------
      # STEIGER
      # --------------------------------------------------------------------------
      
      Steiger_r_exposure = steiger_r_exposure,
      
      Steiger_r_outcome = steiger_r_outcome,
      
      Steiger_p = steiger_p,
      
      Steiger_correct_direction = steiger_direction
      
    )
    
    
    # ============================================================================
    # 9.26 ADD NOMINAL SIGNIFICANCE FLAGS
    # ============================================================================
    
    results <- results %>%
      
      mutate(
        
        IVW_significant_P05 =
          ifelse(
            !is.na(IVW_p) &
              IVW_p < 0.05,
            "Yes",
            "No"
          ),
        
        IVW_significant_P005 =
          ifelse(
            !is.na(IVW_p) &
              IVW_p < 0.005,
            "Yes",
            "No"
          ),
        
        Egger_significant_P05 =
          ifelse(
            !is.na(Egger_p) &
              Egger_p < 0.05,
            "Yes",
            "No"
          ),
        
        WeightedMedian_significant_P05 =
          ifelse(
            !is.na(WeightedMedian_p) &
              WeightedMedian_p < 0.05,
            "Yes",
            "No"
          ),
        
        Weak_instrument_flag =
          ifelse(
            !is.na(Minimum_F) &
              Minimum_F < 10,
            "Yes",
            "No"
          )
        
      )
    
    
    # ============================================================================
    # 9.27 PRINT RESULTS
    # ============================================================================
    
    print(results)
    
    
    # ============================================================================
    # 9.28 STORE RESULT
    # ============================================================================
    
    all_results[[trait]] <- results
    
    
    # ============================================================================
    # 9.29 CLEAN TEMPORARY PLINK FILES
    # ============================================================================
    
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
    
    
    message(
      "Completed: ",
      trait
    )
    
  }
  
  
  # ==============================================================================
  # 10. COMBINE ALL SUMMARY RESULTS
  # ==============================================================================
  
  final_output <- bind_rows(
    all_results
  )
  
  
  # ==============================================================================
  # 11. COMBINE ALL F-STATISTICS
  # ==============================================================================
  
  final_f_statistics <- bind_rows(
    all_f_statistics
  )
  
  
  # ==============================================================================
  # 12. COMBINE ALL HARMONISED DATA
  # ==============================================================================
  
  final_harmonised <- bind_rows(
    all_harmonised_data
  )
  
  
  # ==============================================================================
  # 13. PRINT FINAL RESULTS
  # ==============================================================================
  
  message("")
  message("==============================================================")
  message("FINAL RESULTS")
  message("==============================================================")
  
  
  print(final_output)
  
  
  # ==============================================================================
  # 14. IDENTIFY NOMINALLY SIGNIFICANT IVW RESULTS
  # ==============================================================================
  
  significant_ivw <- final_output %>%
    
    filter(
      !is.na(IVW_p),
      IVW_p < 0.05
    ) %>%
    
    arrange(
      IVW_p
    )
  
  
  message("")
  message("==============================================================")
  message("NOMINALLY SIGNIFICANT IVW RESULTS: P < 0.05")
  message("==============================================================")
  
  
  print(significant_ivw)
  
  
  # ==============================================================================
  # 15. IDENTIFY STRONGER IVW RESULTS
  #
  # P < 0.005 is included as a more stringent threshold.
  # This is NOT a substitute for a formal multiple-testing correction.
  # ==============================================================================
  
  stronger_ivw <- final_output %>%
    
    filter(
      !is.na(IVW_p),
      IVW_p < 0.005
    ) %>%
    
    arrange(
      IVW_p
    )
  
  
  message("")
  message("==============================================================")
  message("IVW RESULTS WITH P < 0.005")
  message("==============================================================")
  
  
  print(stronger_ivw)
  
  
  # ==============================================================================
  # 16. SAVE EXCEL WORKBOOK
  # ==============================================================================
  
  wb <- createWorkbook()
  
  
  # ------------------------------------------------------------------------------
  # Sheet 1: Main Results
  # ------------------------------------------------------------------------------
  
  addWorksheet(
    wb,
    "Main_MR_Results"
  )
  
  writeData(
    wb,
    "Main_MR_Results",
    final_output
  )
  
  
  # ------------------------------------------------------------------------------
  # Sheet 2: Significant IVW
  # ------------------------------------------------------------------------------
  
  addWorksheet(
    wb,
    "Significant_IVW"
  )
  
  writeData(
    wb,
    "Significant_IVW",
    significant_ivw
  )
  
  
  # ------------------------------------------------------------------------------
  # Sheet 3: F Statistics
  # ------------------------------------------------------------------------------
  
  addWorksheet(
    wb,
    "F_Statistics"
  )
  
  writeData(
    wb,
    "F_Statistics",
    final_f_statistics
  )
  
  
  # ------------------------------------------------------------------------------
  # Sheet 4: Harmonised Data
  # ------------------------------------------------------------------------------
  
  addWorksheet(
    wb,
    "Harmonised_Data"
  )
  
  writeData(
    wb,
    "Harmonised_Data",
    final_harmonised
  )
  
  
  # ==============================================================================
  # 17. FORMAT EXCEL
  # ==============================================================================
  
  for (sheet in names(wb)) {
    
    setColWidths(
      wb,
      sheet = sheet,
      cols = 1:100,
      widths = "auto"
    )
    
  }
  
  
  # ==============================================================================
  # 18. SAVE EXCEL
  # ==============================================================================
  
  write.xlsx(
    wb,
    output_file,
    overwrite = TRUE
  )
  
  
  # ==============================================================================
  # 19. FINAL MESSAGE
  # ==============================================================================
  
  message("")
  message("==============================================================")
  message("ANALYSIS COMPLETE")
  message("==============================================================")
  
  message(
    "Outcome: ",
    outcome_name
  )
  
  message(
    "Traits successfully analysed: ",
    length(all_results)
  )
  
  message(
    "Output saved to: ",
    output_file
  )
  
  message("")
  message(
    "Nominally significant IVW results: ",
    nrow(significant_ivw)
  )
  
  message("")
  
  
  # ==============================================================================
  # END OF SCRIPT
  # ==============================================================================
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
  
  out_full <- fread("C:/Users/HP/Documents/msc_project/input/HF_Data/HF_Data/HFpEF_GCST90726619.tsv") %>%
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
    phenotype_col = "HFpEF"
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
  "C:/Users/HP/Documents/msc_project/output/all_LV_traits_IVW_results_HFpEF.xlsx",
  overwrite = TRUE
)

message("All traits processed successfully!")




