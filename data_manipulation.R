# data manipulation

# Create Outputs directory
dir.create("Outputs", showWarnings = FALSE)

# Create directory to store different versions of taxa table.
dir.create('Outputs/Tables/', recursive = T, showWarnings = F)

# Get working files
taxa <- taxa_og
meta <- meta_og

# keep only bacteria in taxa table and set taxa to rownames
taxa <- taxa %>%
  dplyr::filter(str_detect(`@tax`, "Bacteria")) %>%
  tibble::column_to_rownames("@tax")

names(taxa) <- gsub("\\.", "_", names(taxa))

# Check - Make sure all metadata samples are in feature tables and in correct order
stopifnot(all.equal(sort(meta$sample_id), sort(colnames(taxa))))

# Sample-wise filtering of taxa with relative abundance less than 10-5
taxa_relabs <- sweep(taxa, 2, colSums(taxa), "/")
taxa[taxa_relabs < 1e-05] <- 0

# Remove zero abundant taxa
taxa <- taxa[rowSums(taxa) > 0, ]

# Rarefy feature tables to lowest read depth out of all samples
taxa_t <- as.data.frame(t(taxa))
set.seed(725)
taxa_rarefied <- as.data.frame(rrarefy(taxa_t, min(rowSums(taxa_t)))) # 193736


# Filter out samples with less than 1M counts
taxa_clean <- taxa %>%
  dplyr::select(where( ~ sum(.x) > 1E6))

# Filter metadata accordingly and further filter to only keep complete occurrences
meta_clean <- meta %>%
  dplyr::filter(sample_id %in% colnames(taxa_clean)) %>%
  group_by(Animal_Name) %>%
  dplyr::filter(n() == 2) %>%
  ungroup()

# Now keep same samples in taxa table
taxa_clean <- taxa_clean %>%
  dplyr::select(all_of(meta_clean$sample_id))

# Remove zero abundant taxa
taxa_clean <- taxa_clean[rowSums(taxa_clean) > 0, ]

taxa_clean_t <- as.data.frame(t(taxa_clean))
set.seed(725)
taxa_clean_rarefied <- as.data.frame(rrarefy(taxa_clean_t, min(rowSums(taxa_clean_t)))) # 1231130

# Write a copy of the rarefied taxa count table to CSV.
taxa_clean_rare_return <- rownames_to_column(as.data.frame(t(taxa_clean_rarefied)), var = "taxonomy")
write.csv(taxa_clean_rare_return, file = "Outputs/Tables/taxatable_filtered_rarefied.csv", row.names = FALSE)

# Calculate three diversity metrics for rarefied tables
meta_clean$Shannon_OTUs <- diversity(taxa_clean_rarefied, index = "shannon")
meta_clean$Inverse_Simpson_OTUs <- diversity(taxa_clean_rarefied, index = "invsimpson")
meta_clean$Observed_OTUs <- rowSums(taxa_clean_rarefied != 0)

# Convert rarefied tables to relative abundances, and re-transpose for vegdist for beta diversity plotting later
taxa_clean_rarefied_relative <- as.data.frame(apply(taxa_clean_rarefied, 1, function(x) x/sum(x)))
taxa_clean_rarefied_relative_t <- as.data.frame(t(taxa_clean_rarefied_relative))
