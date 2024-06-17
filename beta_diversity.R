##### beta diversity #####
library(vroom)
# need to specify the beta tables that will be analyzed 
beta_tables <- c("taxa", "koe", "vfs", "amr", "caz")

# Create outputs directory
output_dir_beta <- "Outputs/beta_diversity"
set.seed(513)

#long tcam and permanova

for (beta_table in beta_tables) {
  output_dir_beta_long_permanova <- file.path(output_dir_beta,"longitudinal/permanovaFL",beta_table)
  output_dir_beta_long_tcam <- file.path(output_dir_beta,"longitudinal/tcam",beta_table)
  dir.create(output_dir_beta_long_tcam, recursive = T, showWarnings = F)
  dir.create(output_dir_beta_long_permanova, recursive = T, showWarnings = F)
  
}


# Create a dataframe to store statistical results from permanovas.
permanovafl_data <- data.frame(data_type = character(), term = character(), `F` = numeric(), R2 = numeric(), P = numeric())


for (beta_table in beta_tables){
  if(beta_table == "taxa"){
    meta_tmp <- meta_taxa
    plot_label <- "Taxa"
  }else if(beta_table == "koe"){
    meta_tmp <- meta_koe
    plot_label <- "KEGG Enzymes"
  }else if(beta_table == "vfs"){
    meta_tmp <- meta_vfs
    plot_label <- "Virulence factors"
  }else if(beta_table == "amr"){
    meta_tmp <- meta_amr
    plot_label <- "Antimicrobial resistance genes"
  }else if(beta_table == "caz"){
    meta_tmp <- meta_caz
    plot_label <- "CAZymes"
  }
  tmp_DM <- as.data.frame(get.clr(get(beta_table)))
  DM_ <- as.matrix(vegdist(tmp_DM, 'euclidean'))
  DM <- as.data.frame(DM_)
  write.csv(DM, paste0("Outputs/beta_diversity/longitudinal/permanovaFL/",beta_table,'/',(beta_table), '_permanova_dist.csv'))
  
  # order metadata to be in the same order as data.
  meta_tmp$Sample_ID <- sub('_', '.', meta_tmp$Sample_ID)
  
  meta_tmp <- meta_tmp[match(row.names(DM), meta_tmp$Sample_ID),]
  
  # make sure sample names are in the correct order.
  stopifnot(all.equal(row.names(DM), sub('_', '.', meta_tmp$Sample_ID)))
  set.seed(513)
  permanovaFL_out <- permanovaFL(DM ~ Timepoint + Treatment:Timepoint,
                                 data = as.data.frame(meta_tmp),
                                 cluster.id = Subject_ID,
                                 perm.within.type = "free",
                                 perm.between.type = "none",
                                 n.perm.max = nperm)
  
  # Get permanovaFL output in plotting format.
  permanovaFL_F_interaction <- toString(round(permanovaFL_out$F.statistics[2],4))
  permanovaFL_R2_interaction <- toString(round(permanovaFL_out$R.squared[2],4))
  permanovaFL_P_interaction <- toString(p_format(permanovaFL_out$p.permanova[2],4, accuracy = 1e-04))
  permanovaFL_F_timepoint <- toString(round(permanovaFL_out$F.statistics[1],4))
  permanovaFL_R2_timepoint <- toString(round(permanovaFL_out$R.squared[1],4))
  permanovaFL_P_timepoint <- toString(round(permanovaFL_out$p.permanova[1],4))
  
  permanovaFL_interaction_output <- paste("permanovaFL Treatment Group x Time Point: \nF = ", permanovaFL_F_interaction, "; R2 = ", permanovaFL_R2_interaction, "; p = ", permanovaFL_P_interaction, sep = "")
  permanovaFL_timepoint_output <- paste("permanovaFL Timepoint: F = ", permanovaFL_F_timepoint, "; R2 = ", permanovaFL_R2_timepoint, "; p = ", permanovaFL_P_timepoint, sep = "")
  # Add permanovaFL output to statistics data.frame.
  permanovafl_data_tmp <- data.frame(data_type = beta_table, term = c('Time_Point', 'Treatment:Time_Point'), `F` = permanovaFL_out$F.statistics, R2 = permanovaFL_out$R.squared, P = permanovaFL_out$p.permanova)
  
  # grow statistics dataframe.
  permanovafl_data <- rbind(permanovafl_data, permanovafl_data_tmp)
  # Make a PCoA object.
  pcoa_object <- ape::pcoa(DM)
  
  pc_vectors <- pcoa_object$vectors
  pc_1_and_2 <- pc_vectors[,1:2]
  pc_1_and_2 <- data.frame(pc_1_and_2)
  colnames(pc_1_and_2) <- c("PC1", "PC2")
  # Grab % explained
  pc_vals <- pcoa_object$values$Relative_eig
  pc_val_1 <- round(pc_vals[1]*100, digits = 2)
  pc_val_2 <- round(pc_vals[2]*100, digits = 2)
  # Add metadata.
  pc_1_and_2$Timepoint <- meta_tmp$Timepoint
  pc_1_and_2$Treatment <- meta_tmp$Treatment
  pc_1_and_2$Subject_ID <- meta_tmp$Subject_ID
  # make treatment_group_timepoint combined variable for color scale plotting
  pc_1_and_2$Treatment_Group_Timepoint <- paste(pc_1_and_2$Treatment, pc_1_and_2$Timepoint, sep = ' ')
  pc_1_and_2$Diet_Group_Timepoint <- factor(pc_1_and_2$Treatment_Group_Timepoint, levels = c("HM 1","SF 1","CF 1","SF 2","CF 2","HM 2"))
  
  # Make PCoA plots.
  # geom_path(aes(group = .data[['Subject_ID']]), color = 'grey', alpha = 0.5, arrow= arrow(type = "open", length = unit(0.2, "cm")))
  # Create the treatment:timepoint plot.
  # Color scheme found here: https://meyerweb.com/eric/tools/color-blend/#FF0000:FFCC33:2:hex
  tmp_plot  <- ggplot(data = pc_1_and_2, aes(x = .data[["PC1"]], 
                                             y = .data[["PC2"]], color = .data[["Treatment"]], shape = .data[["Timepoint"]])) +
    geom_point(size = 3, alpha = 0.7)  + theme(legend.position = "right",  
                                               legend.text = element_text(size = 14),
                                               legend.title = element_blank(), 
                                               plot.caption = element_text(hjust = 0)) +
    scale_color_manual(values = abbott_colors) +
    xlab(paste("PC1 (", pc_val_1, "%)", sep = "")) +
    ylab(paste("PC2 (", pc_val_2, "%)", sep = "")) +
    # labs(caption = paste0(permanovaFL_timepoint_output, '\n', permanovaFL_interaction_output)) +
    labs(caption = permanovaFL_interaction_output) +
    ggtitle(paste0('Treatment groups over time')) +
    labs(subtitle = paste0(plot_label,': Aitchiston distance'))
  
  
  # # Save plots.
  pdf(paste0('Outputs/beta_diversity/longitudinal/permanovaFL/', beta_table, "/PCA_treatment_Timepoint_Interaction.pdf"), width = 5.7, height = 5)
  print(tmp_plot)
  dev.off()
  
}

# Save permanovaFL data as a .csv file.
write.csv(permanovafl_data, 'Outputs/beta_diversity/longitudinal/permanovaFL/permanovaFL_results.csv', row.names = FALSE)

###################################################################################################################FALSE#######################
#	Longitudinal Approach -- TCAM

#output_dir_beta_long_tcam

for (beta_table in beta_tables){
  if(beta_table == "taxa"){
    meta_tmp <- meta_taxa
    plot_label <- "Taxa"
  }else if(beta_table == "koe"){
    meta_tmp <- meta_koe
    plot_label <- "KEGG Enzymes"
  }else if(beta_table == "vfs"){
    meta_tmp <- meta_vfs
    plot_label <- "Virulence factors"
  }else if(beta_table == "amr"){
    meta_tmp <- meta_amr
    plot_label <- "Antimicrobial resistance genes"
  }else if(beta_table == "caz"){
    meta_tmp <- meta_caz
    plot_label <- "CAZymes"
  }
  output_dir_beta_long_tcam <- file.path(output_dir_beta,"longitudinal/tcam",beta_table)
  # subset to remove 0 abundance taxa.
  b_table <- get(beta_table)
  min(rowSums(b_table))
  b_table <- b_table[rowSums(b_table) > 0,]
  tmp_DM <- as.data.frame(get.clr(b_table))
  meta_tmp <- meta_tmp %>% dplyr::select(c(Sample_ID, Timepoint,Treatment,Subject_ID))
  rownames(tmp_DM) <- gsub("\\.", "_", rownames(tmp_DM))
  # join metadata and transformed taxa table.
  #clr_meta <- merge(meta_tmp, tmp_DM, by.x = 'Sample_Id', by.y = 'row.names')
  feature_table_clr_meta <- merge(meta_tmp, tmp_DM, by.x = 'Sample_ID', by.y = 'row.names')
  feature_table_clr_meta <- feature_table_clr_meta %>%
    dplyr::arrange(Subject_ID, Timepoint) %>%
    dplyr::group_by(Timepoint) 
  
  # move sample id to row.names
  #need to ensure every 
  
  feature_table_clr_meta <- feature_table_clr_meta %>%
    group_by(Subject_ID) %>%
    filter(all(c(1, 2) %in% Timepoint)) %>%
    ungroup()
  
  feature_table_clr_meta <- feature_table_clr_meta %>% column_to_rownames('Sample_ID')
  feature_table_clr_meta$Timepoint <- factor(feature_table_clr_meta$Timepoint, levels = unique(feature_table_clr_meta$Timepoint))
  
  #### Need to do calculations for deviation from baseline
  feature_table_clr_meta_melted <- 
    feature_table_clr_meta %>%reshape2::melt(id.vars = c('Subject_ID', 'Timepoint',"Treatment"))
  
  #Have two time points 2 and 1, assuming 1 is the baseline
  feature_table_clr_meta_melted$Timepoint <- as.numeric(feature_table_clr_meta_melted$Timepoint)
  feature_table_clr_meta_melted$value <- as.numeric(feature_table_clr_meta_melted$value)
  
  feature_table_clr_meta_melted_dfb <- feature_table_clr_meta_melted %>%
    dplyr::group_by(Subject_ID, variable) %>%
    dplyr::mutate(DFB = value - value[Timepoint == 1]) %>%
    dplyr::filter(Timepoint != 1)
  
  # Now convert to wide format using DFB values for the TCAM algorithm to work with.
  feature_table_for_tcam_dfb <- feature_table_clr_meta_melted_dfb %>%
    dplyr::select(-value) %>%
    tidyr::pivot_wider(id_cols = c('Subject_ID', 'Timepoint'), values_from = 'DFB', names_from = 'variable') %>%
    # arrange by subject_id and week
    dplyr::arrange(Subject_ID, Timepoint)
  
  vroom_write(feature_table_for_tcam_dfb,paste0(output_dir_beta_long_tcam,'/table_dfb_meta.tsv'))
  #need to run code that runs TCAM.py
  
}

system('source code/run_TCAM.sh')

# TCAM.R part 2

for (beta_table in beta_tables){
  if(beta_table == "taxa"){
    meta_tmp <- meta_taxa
    plot_label <- "Taxa"
  }else if(beta_table == "koe"){
    meta_tmp <- meta_koe
    plot_label <- "KEGG Enzymes"
  }else if(beta_table == "vfs"){
    meta_tmp <- meta_vfs
    plot_label <- "Virulence factors"
  }else if(beta_table == "amr"){
    meta_tmp <- meta_amr
    plot_label <- "Antimicrobial resistance genes"
  }else if(beta_table == "caz"){
    meta_tmp <- meta_caz
    plot_label <- "CAZymes"
  }
  
  # Read in output from python script.
  # This should have Subjects (subject_id) as the rows and factors as columns.
  # From this I can calculate Euclidean distances between Subjects.
  output_dir_beta_long_tcam <- file.path(output_dir_beta,"longitudinal/tcam",beta_table)
  tca_transform <- vroom(paste0(output_dir_beta_long_tcam, '/',(beta_table), '_tcam_transform.csv'))
  tca_transform <- column_to_rownames(tca_transform, '...1')
  # Make a distance matrix based on factors
  tca_transform_dist <- as.matrix(dist(tca_transform, method = 'euclidean'))
  # Save the distance matrix.
  write.csv(tca_transform_dist, paste0(output_dir_beta_long_tcam,'/', toupper(beta_table), '_tcam_dist.csv'))
  
  # subset metadata to the subjects in the TCAM matrix, and just use the metadata for one datapoint per subject.
  meta_subjects <- meta_tmp %>% 
    dplyr::filter(Subject_ID %in% row.names(tca_transform_dist)) %>%
    dplyr::group_by(Subject_ID) %>%
    dplyr::slice_head(n = 1)
  
  # Make sure metadata samples are in the same order.
  stopifnot(all.equal(as.character(meta_subjects$Subject_ID), row.names(tca_transform_dist)))
  permanova_out <- adonis2(tca_transform_dist ~ Treatment, data = meta_subjects)
  # Save PERMANOVA results.
  utils::capture.output(permanova_out, file = paste0(output_dir_beta_long_tcam, '/', toupper(beta_table), '_tcam_permanova.txt'))
  
  # Make a TCA plot.
  tca_object <- tca_transform
  
  tca_1_and_2 <- tca_object[,1:2]
  tca_1_and_2 <- data.frame(tca_1_and_2)
  colnames(tca_1_and_2) <- c("F1", "F2")
  # Grab % explained
  tca_percent_explained <- colnames(vroom(paste0(output_dir_beta_long_tcam, '/', toupper(beta_table), '_tcam_variance_explained.csv')))
  F_val_1 <- round(as.numeric(tca_percent_explained[1]), digits = 2)
  F_val_2 <- round(as.numeric(tca_percent_explained[2]), digits = 2)
  
  # Add metadata.
  tca_1_and_2$Treatment <- meta_subjects$Treatment
  tca_1_and_2$Subject_ID <- meta_subjects$Subject_ID
  # Create string for display on plot.
  permanova_output_string <- paste("PERMANOVA\nTreatment Group: F = ", round(permanova_out$F[1], 3), "; R2 = ", round(permanova_out$R2[1], 3), "; p = ", p_format(permanova_out$`Pr(>F)`[1], accuracy = 1e-04), sep = "")
  
  tca_1_and_2$Treatment <- as.character(tca_1_and_2$Treatment)
  
  tcam_plot <- ggplot(data = tca_1_and_2, aes(x = .data[["F1"]], 
                                              y = .data[["F2"]],
                                              color = .data[["Treatment"]])) +
    geom_point(size = 3, alpha = 0.7) +
    theme(legend.position = "right",  
          legend.text = element_text(size = 14),
          legend.title = element_blank(), 
          plot.caption = element_text(hjust = 0)) +
    scale_color_manual(values = abbott_colors) + 
    xlab(paste("F1 (", F_val_1, "%)", sep = "")) +
    ylab(paste("F2 (", F_val_2, "%)", sep = "")) +
    labs(caption = permanova_output_string) +
    ggtitle(paste0('Treatment group over time')) +
    labs(subtitle = toupper(plot_label))
  pdf(paste0(output_dir_beta_long_tcam,'/', toupper(beta_table), '_TCAM_Treatment_Group.pdf'), width = 5, height = 4.25)
  print(tcam_plot)
  dev.off()
  
}
