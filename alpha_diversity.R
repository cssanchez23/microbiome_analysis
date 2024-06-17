# Alpha diversity

# Make outputs directory
alpha_out_dir <- file.path("Outputs/Alpha_Diversity")
dir.create(path = alpha_out_dir, recursive = TRUE, showWarnings = FALSE)

# Generate alpha outputs table and write to include in deliberables
meta_alpha <- data.frame(Sample_ID = meta_clean$sample_id,
                         Animal_name = meta_clean$Animal_Name,
                         Timepoint = as.factor(meta_clean$Timepoint),
                         Treatment = as.factor(meta_clean$Treatment_Group),
                         Shannon_OTUs = meta_clean$Shannon_OTUs,
                         Inverse_Simpson_OTUs = meta_clean$Inverse_Simpson_OTUs,
                         Observed_OTUs = meta_clean$Observed_OTUs)

# Write out as a deliverable
write.csv(meta_alpha, "Outputs/Tables/alpha_diversity_values.csv", row.names = FALSE)

# Create empty table for storing results of statistical testing
lmer_results <- data.frame(Feature = character(),
                           Alpha_Measure = character(),
                           Term_or_Contrast = character(),
                           Estimate = numeric(),
                           std_error = numeric(),
                           df = numeric(),
                           t_value = numeric(),
                           p_value = numeric())

alpha_measures <- c("Shannon_OTUs", "Inverse_Simpson_OTUs", "Observed_OTUs")

for (alpha_measure in alpha_measures) {
  
  alpha_temp <- meta_alpha
  
  feature <- "OTUs"
  
  alpha_name <- ifelse(str_detect(alpha_measure, "Observed"), "Observed_Features",
                       ifelse(str_detect(alpha_measure, "Shannon"), "Shannon", "Inverse_Simpson"))
  
  #run model
  mod <- lmerTest::lmer(get(alpha_measure) ~ Treatment*Timepoint + (1|Animal_name), data = alpha_temp, REML = FALSE)
 
  #extract coefficients
  mod_coef <- as.data.frame(round(coef(summary(mod)), 4))
  
  # Get coefficients
  coef_temp <- as.data.frame(coef(summary(mod))) %>%
    dplyr::filter(str_detect(rownames(.), "Treatment"))
  
  # Populate results table and merge
  lmer_results_temp <- data.frame(Feature = feature,
                                  Alpha_Measure = alpha_name,
                                  Term_or_Contrast = rownames(coef_temp),
                                  Estimate = coef_temp$Estimate,
                                  std_error = coef_temp$`Std. Error`,
                                  df = coef_temp$df,
                                  t_value = coef_temp$`t value`,
                                  p_value = coef_temp$`Pr(>|t|)`)
  
  lmer_results <- rbind(lmer_results, lmer_results_temp)
  
  # Now calculate the contrasts along with p-values
  
  #contrasts for treatments + timepoint
  mod_con <- as.data.frame(lmerTest::difflsmeans(mod, test.effs = "Treatment")) %>%
    tibble::rownames_to_column("contrast")

  mod_con_out <- mod_con
  
  # first separate the contrast terms from the effect size table
  mod_con$contrast1 <- str_split(mod_con$contrast, '-', simplify = TRUE)[,1]
  mod_con$contrast2 <- str_split(mod_con$contrast, '-', simplify = TRUE)[,2]
  
  contrast_results_out_temp <- data.frame(Feature = feature,
                                          Alpha_Measure = alpha_name,
                                          Term_or_Contrast = mod_con_out$contrast,
                                          Estimate = mod_con_out$Estimate,
                                          std_error = mod_con_out$`Std. Error`,
                                          df = mod_con_out$df,
                                          t_value = mod_con_out$`t value`,
                                          p_value = mod_con_out$`Pr(>|t|)`)
  
  lmer_results <- rbind(lmer_results, contrast_results_out_temp)

}

# Keep only relevant terms
relevant_terms <- c("Treatment5650:TimepointFinal", 
                    "Treatment7687:TimepointFinal",
                    "Treatment2927 - Treatment5650", 
                    "Treatment2927 - Treatment7687", 
                    "Treatment5650 - Treatment7687",
                    "Treatment2927:TimepointBaseline - Treatment5650:TimepointBaseline",
                    "Treatment2927:TimepointBaseline - Treatment7687:TimepointBaseline",
                    "Treatment5650:TimepointBaseline - Treatment7687:TimepointBaseline",
                    "Treatment2927:TimepointFinal - Treatment5650:TimepointFinal",
                    "Treatment2927:TimepointFinal - Treatment7687:TimepointFinal",
                    "Treatment5650:TimepointFinal - Treatment7687:TimepointFinal")

# Adjust p-values
lmer_results_out <- lmer_results %>%
  dplyr::filter(!str_detect(Feature, "Other")) %>%
  dplyr::filter(Term_or_Contrast %in% relevant_terms) %>%
  as.data.frame()

write.csv(lmer_results_out, paste0(
  "Outputs/Alpha_Diversity/Linear_mixed_effects_model_results.csv"
), row.names = FALSE)


