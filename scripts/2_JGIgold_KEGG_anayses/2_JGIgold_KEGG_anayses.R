#!/usr/bin/env Rscript
# JGI/GOLD KEGG Analysis: Bacterial Genome Size Constraints
#
# This script implements comprehensive KEGG annotation and feature extraction:
# 1. Compute genome metrics (coding density, gene counts, etc.)
# 2. Extract amino acid composition and nitrogen burden
# 3. Identify transcription factors and mobile elements
# 4. Prepare for KEGG annotation (KofamScan or KEGG API)
# 5. Map KOs to modules and pathways
# 6. Create comprehensive feature matrix for downstream modeling
#
# OUTPUT: All results saved to results/2_JGIgold_KEGG_anayses_out/

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Biostrings)
  library(ggplot2)
})

# Set paths
BASE_DIR <- "/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
RESULTS_DIR <- file.path(BASE_DIR, "results", "2_JGIgold_KEGG_anayses_out")
ASSEMBLIES_DIR <- file.path(BASE_DIR, "data", "ncbi", "assemblies")

# Create output directory
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Input data paths
NCBI_METADATA <- file.path(BASE_DIR, "data", "ncbi", "metadata", "assembly_data_report_extracted.tsv")
NCBI_BASIC <- file.path(BASE_DIR, "data", "ncbi", "metadata", "metadata_table.tsv")
SELECTED_ENV <- file.path(BASE_DIR, "results", "1_exploratory_analyses_out", "08_selected_environments_summary.tsv")

cat("================================================================================\n")
cat("JGI/GOLD KEGG ANALYSIS: Bacterial Genome Size Constraints\n")
cat("================================================================================\n")
cat("\nOutput directory:", RESULTS_DIR, "\n")
cat("================================================================================\n\n")

# ============================================================================
# 1. Load and Prepare Data
# ============================================================================
cat("[1] Loading metadata and preparing data...\n")

ncbi_detailed <- fread(NCBI_METADATA, sep = "\t", data.table = FALSE)
ncbi_basic <- fread(NCBI_BASIC, sep = "\t", data.table = FALSE)

cat("  ✓ NCBI detailed metadata:", nrow(ncbi_detailed), "genomes\n")
cat("  ✓ NCBI basic metadata:", nrow(ncbi_basic), "genomes\n")

# Load selected environments if available
selected_envs <- NULL
if (file.exists(SELECTED_ENV)) {
  selected_envs <- fread(SELECTED_ENV, sep = "\t", data.table = FALSE)
  cat("  ✓ Selected environments:", nrow(selected_envs), "genomes\n")
}

# Get list of accessions to process
if (!is.null(selected_envs) && "accession" %in% colnames(selected_envs)) {
  accessions <- selected_envs$accession
  cat("  Processing", length(accessions), "genomes from selected environments\n")
} else {
  accessions <- ncbi_detailed$accession[!is.na(ncbi_detailed$accession)]
  cat("  Processing all", length(accessions), "genomes\n")
}

# Limit to first 50 for testing (remove this limit for full run)
# accessions <- accessions[1:50]
cat("  Processing", length(accessions), "genomes\n")

# ============================================================================
# 2. Compute Genome Metrics
# ============================================================================
cat("\n[2] Computing genome metrics...\n")
cat("  Computing: coding density, gene counts, GC content\n")

compute_genome_metrics <- function(acc) {
  acc_dir <- file.path(ASSEMBLIES_DIR, acc)
  
  if (!dir.exists(acc_dir)) {
    return(NULL)
  }
  
  # Get metadata for this accession
  meta_idx <- which(ncbi_detailed$accession == acc)
  if (length(meta_idx) == 0) {
    return(NULL)
  }
  meta_row <- ncbi_detailed[meta_idx[1], ]
  
  # Helper function for safe extraction
  safe_extract <- function(x, default = 0) {
    if (is.null(x) || is.na(x) || x == "") return(default)
    as.numeric(x)
  }
  
  metrics <- data.frame(
    accession = acc,
    genome_size_bp = safe_extract(meta_row$stats_totalSequenceLength),
    gc_percent = safe_extract(meta_row$stats_gcPercent),
    genes_total = safe_extract(meta_row$genes_total),
    genes_proteinCoding = safe_extract(meta_row$genes_proteinCoding),
    checkm_completeness = safe_extract(meta_row$checkm_completeness),
    checkm_contamination = safe_extract(meta_row$checkm_contamination),
    stringsAsFactors = FALSE
  )
  
  # Compute coding density from CDS file
  cds_file <- file.path(acc_dir, "cds_from_genomic.fna")
  if (file.exists(cds_file)) {
    tryCatch({
      cds_seqs <- readDNAStringSet(cds_file)
      total_cds_length <- sum(width(cds_seqs))
      if (metrics$genome_size_bp > 0) {
        metrics$coding_density <- (total_cds_length / metrics$genome_size_bp) * 100
      } else {
        metrics$coding_density <- NA
      }
    }, error = function(e) {
      metrics$coding_density <- NA
    })
  } else {
    metrics$coding_density <- NA
  }
  
  return(metrics)
}

# Process genomes in batches
genome_metrics_list <- list()
batch_size <- 50  # Larger batches for efficiency

for (i in seq(1, length(accessions), by = batch_size)) {
  batch <- accessions[i:min(i + batch_size - 1, length(accessions))]
  batch_results <- lapply(batch, compute_genome_metrics)
  batch_results <- batch_results[!sapply(batch_results, is.null)]
  genome_metrics_list <- c(genome_metrics_list, batch_results)
  
  if (i %% 100 == 1 || i + batch_size > length(accessions)) {
    cat("    Processed", length(genome_metrics_list), "genomes...\n")
  }
}

genome_metrics_df <- bind_rows(genome_metrics_list)
cat("  ✓ Computed metrics for", nrow(genome_metrics_df), "genomes\n")

# Save genome metrics
fwrite(genome_metrics_df, file.path(RESULTS_DIR, "01_genome_metrics.tsv"), sep = "\t")
cat("  ✓ Saved:", file.path(RESULTS_DIR, "01_genome_metrics.tsv"), "\n")

# ============================================================================
# 3. Amino Acid Composition and Nitrogen Burden
# ============================================================================
cat("\n[3] Computing amino acid composition and nitrogen burden...\n")
cat("  Rationale: Nitrogen burden indicates nutrient limitation signatures\n")

# Amino acid nitrogen content (atoms per residue)
AA_N_CONTENT <- c(
  A = 1, R = 4, N = 2, D = 1, C = 1, Q = 2, E = 1, G = 1,
  H = 3, I = 1, L = 1, K = 2, M = 1, F = 1, P = 1, S = 1,
  T = 1, W = 2, Y = 1, V = 1, X = 0, `*` = 0
)

compute_aa_composition <- function(acc) {
  acc_dir <- file.path(ASSEMBLIES_DIR, acc)
  protein_file <- file.path(acc_dir, "protein.faa")
  
  if (!file.exists(protein_file)) {
    return(NULL)
  }
  
  tryCatch({
    protein_seqs <- readAAStringSet(protein_file)
    
    if (length(protein_seqs) == 0) {
      return(NULL)
    }
    
    # Concatenate all protein sequences
    all_aa <- paste(as.character(protein_seqs), collapse = "")
    all_aa <- toupper(all_aa)
    total_aa <- nchar(all_aa)
    
    if (total_aa == 0) {
      return(NULL)
    }
    
    # Count amino acids
    aa_counts <- table(strsplit(all_aa, "")[[1]])
    
    # Calculate nitrogen burden (average N atoms per residue)
    total_n_atoms <- sum(AA_N_CONTENT[names(aa_counts)] * aa_counts, na.rm = TRUE)
    n_burden <- total_n_atoms / total_aa
    
    # Calculate amino acid frequencies
    aa_freq <- (aa_counts / total_aa) * 100
    aa_freq_df <- as.data.frame(t(aa_freq))
    colnames(aa_freq_df) <- paste0("aa_", colnames(aa_freq_df))
    
    aa_data <- data.frame(
      accession = acc,
      amino_n_burden = n_burden,
      total_aa_residues = total_aa,
      stringsAsFactors = FALSE
    )
    aa_data <- cbind(aa_data, aa_freq_df)
    
    return(aa_data)
  }, error = function(e) {
    return(NULL)
  })
}

# Process amino acid composition
aa_list <- list()
for (i in seq(1, length(accessions), by = batch_size)) {
  batch <- accessions[i:min(i + batch_size - 1, length(accessions))]
  batch_results <- lapply(batch, compute_aa_composition)
  batch_results <- batch_results[!sapply(batch_results, is.null)]
  aa_list <- c(aa_list, batch_results)
  
  if (i %% 100 == 1 || i + batch_size > length(accessions)) {
    cat("    Processed", length(aa_list), "genomes...\n")
  }
}

# Combine amino acid composition data
if (length(aa_list) > 0) {
  aa_df <- bind_rows(aa_list)
  cat("  ✓ Computed amino acid composition for", nrow(aa_df), "genomes\n")
  
  # Save amino acid composition
  fwrite(aa_df, file.path(RESULTS_DIR, "02_amino_acid_composition.tsv"), sep = "\t")
  cat("  ✓ Saved:", file.path(RESULTS_DIR, "02_amino_acid_composition.tsv"), "\n")
} else {
  aa_df <- data.frame(accession = character(), amino_n_burden = numeric(), stringsAsFactors = FALSE)
  cat("  ⚠ No amino acid composition data computed\n")
}

# ============================================================================
# 4. Transcription Factor and Mobile Element Counts
# ============================================================================
cat("\n[4] Identifying transcription factors and mobile elements...\n")
cat("  Rationale: TF count indicates regulatory complexity\n")
cat("  Mobile elements indicate horizontal gene transfer\n")

# Keywords for transcription factors and mobile elements
TF_KEYWORDS <- c(
  "transcription factor", "transcriptional regulator", "helix-turn-helix",
  "sigma factor", "response regulator", "lysr", "arac", "gntr", "iclr",
  "marr", "tetr", "merr", "luxr", "fis", "hns", "ihf", "fru", "crp"
)

MOBILE_KEYWORDS <- c(
  "transposase", "transposon", "integrase", "recombinase", "conjugative",
  "mobilization", "plasmid", "phage", "prophage", "insertion sequence",
  "is element", "tnp", "ins", "tra", "mob"
)

count_tf_mobile <- function(acc) {
  acc_dir <- file.path(ASSEMBLIES_DIR, acc)
  gff_files <- list.files(acc_dir, pattern = "\\.gff$", full.names = TRUE)
  
  if (length(gff_files) == 0) {
    return(NULL)
  }
  
  gff_file <- gff_files[1]
  
  tryCatch({
    gff_lines <- readLines(gff_file)
    gff_lines <- gff_lines[!grepl("^#", gff_lines)]
    
    tf_count <- 0
    mobile_count <- 0
    
    for (line in gff_lines) {
      line_lower <- tolower(line)
      
      # Check for transcription factors
      if (any(sapply(TF_KEYWORDS, function(kw) grepl(kw, line_lower, fixed = TRUE)))) {
        tf_count <- tf_count + 1
      }
      
      # Check for mobile elements
      if (any(sapply(MOBILE_KEYWORDS, function(kw) grepl(kw, line_lower, fixed = TRUE)))) {
        mobile_count <- mobile_count + 1
      }
    }
    
    return(data.frame(
      accession = acc,
      tf_count = tf_count,
      mobile_element_count = mobile_count,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    return(NULL)
  })
}

# Process TF and mobile element counts
tf_mobile_list <- list()
for (i in seq(1, length(accessions), by = batch_size)) {
  batch <- accessions[i:min(i + batch_size - 1, length(accessions))]
  batch_results <- lapply(batch, count_tf_mobile)
  batch_results <- batch_results[!sapply(batch_results, is.null)]
  tf_mobile_list <- c(tf_mobile_list, batch_results)
  
  if (i %% 100 == 1 || i + batch_size > length(accessions)) {
    cat("    Processed", length(tf_mobile_list), "genomes...\n")
  }
}

if (length(tf_mobile_list) > 0) {
  tf_mobile_df <- bind_rows(tf_mobile_list)
  cat("  ✓ Identified TFs and mobile elements for", nrow(tf_mobile_df), "genomes\n")
  
  # Save TF and mobile element counts
  fwrite(tf_mobile_df, file.path(RESULTS_DIR, "03_tf_mobile_elements.tsv"), sep = "\t")
  cat("  ✓ Saved:", file.path(RESULTS_DIR, "03_tf_mobile_elements.tsv"), "\n")
} else {
  tf_mobile_df <- data.frame(accession = character(), tf_count = numeric(), 
                            mobile_element_count = numeric(), stringsAsFactors = FALSE)
  cat("  ⚠ No TF/mobile element data computed\n")
}

# ============================================================================
# 5. KEGG Annotation Preparation
# ============================================================================
cat("\n[5] Preparing for KEGG annotation...\n")
cat("  Note: KEGG annotation can be done using:\n")
cat("    - KofamScan (command-line tool) - recommended for large-scale annotation\n")
cat("    - KEGG API (R package KEGGREST) - for querying KEGG database\n")
cat("    - KEGG Mapper for module mapping\n")

# Check if KEGGREST is available
keggrestr_available <- requireNamespace("KEGGREST", quietly = TRUE)

if (keggrestr_available) {
  library(KEGGREST)
  cat("  ✓ KEGGREST package available\n")
  cat("    Note: KEGGREST can query KEGG database but KofamScan is needed for annotation\n")
} else {
  cat("  ⚠ KEGGREST not installed - install with: BiocManager::install('KEGGREST')\n")
}

# For now, create placeholder KEGG results
# In production, this would be populated from KofamScan output
cat("  Creating placeholder KEGG results (to be updated with KofamScan output)\n")

# Use processed accessions for KEGG results
processed_accs <- unique(genome_metrics_df$accession)
kegg_df <- data.frame(
  accession = processed_accs,
  ko_count = NA_real_,
  module_count = NA_real_,
  ko_per_mb = NA_real_,
  module_per_mb = NA_real_,
  stringsAsFactors = FALSE
)

# Save KEGG results placeholder
fwrite(kegg_df, file.path(RESULTS_DIR, "04_kegg_annotation.tsv"), sep = "\t")
cat("  ✓ Saved placeholder:", file.path(RESULTS_DIR, "04_kegg_annotation.tsv"), "\n")
cat("    (Update with actual KEGG annotation results from KofamScan)\n")

# ============================================================================
# 6. Merge All Features into Comprehensive Matrix
# ============================================================================
cat("\n[6] Creating comprehensive feature matrix...\n")

# Start with processed accessions and get their metadata
processed_accessions <- unique(genome_metrics_df$accession)

# Get metadata for processed genomes only
feature_matrix <- ncbi_basic %>%
  filter(`Assembly Accession` %in% processed_accessions) %>%
  select(`Assembly Accession`, `Organism Name`, `Organism Taxonomic ID`,
         `Assembly Release Date`, `Assembly Level`) %>%
  rename(accession = `Assembly Accession`,
         organism_name = `Organism Name`,
         organism_taxId = `Organism Taxonomic ID`,
         release_date = `Assembly Release Date`,
         assembly_level = `Assembly Level`)

# Merge genome metrics
feature_matrix <- feature_matrix %>%
  inner_join(genome_metrics_df, by = "accession")

# Merge amino acid composition
if (nrow(aa_df) > 0 && "amino_n_burden" %in% colnames(aa_df)) {
  feature_matrix <- feature_matrix %>%
    left_join(aa_df %>% select(accession, amino_n_burden), by = "accession")
} else {
  feature_matrix$amino_n_burden <- NA
}

# Merge TF and mobile element counts
if (nrow(tf_mobile_df) > 0) {
  feature_matrix <- feature_matrix %>%
    left_join(tf_mobile_df, by = "accession")
} else {
  feature_matrix$tf_count <- NA
  feature_matrix$mobile_element_count <- NA
}

# Merge KEGG results
feature_matrix <- feature_matrix %>%
  left_join(kegg_df, by = "accession")

# Add environment information if available
if (!is.null(selected_envs) && "ORGANISM ECOSYSTEM CATEGORY" %in% colnames(selected_envs)) {
  env_map <- selected_envs %>%
    select(accession, `ORGANISM ECOSYSTEM CATEGORY`) %>%
    distinct()
  feature_matrix <- feature_matrix %>%
    left_join(env_map %>% rename(environment = `ORGANISM ECOSYSTEM CATEGORY`), 
              by = "accession")
}

# Compute derived metrics
feature_matrix <- feature_matrix %>%
  mutate(
    genome_size_mb = genome_size_bp / 1e6,
    ko_per_mb = ko_count / genome_size_mb,
    module_per_mb = module_count / genome_size_mb
  )

# Save comprehensive feature matrix
fwrite(feature_matrix, file.path(RESULTS_DIR, "05_genome_feature_matrix.tsv"), sep = "\t")
cat("  ✓ Saved comprehensive feature matrix:", 
    file.path(RESULTS_DIR, "05_genome_feature_matrix.tsv"), "\n")
cat("    ", nrow(feature_matrix), "genomes,", ncol(feature_matrix), "features\n")

# ============================================================================
# 7. Summary Statistics
# ============================================================================
cat("\n[7] Computing summary statistics...\n")

summary_stats <- data.frame(
  total_genomes = nrow(feature_matrix),
  genomes_with_metrics = sum(!is.na(feature_matrix$genome_size_bp)),
  genomes_with_aa = sum(!is.na(feature_matrix$amino_n_burden)),
  genomes_with_tf = sum(!is.na(feature_matrix$tf_count)),
  genomes_with_mobile = sum(!is.na(feature_matrix$mobile_element_count)),
  genomes_with_kegg = sum(!is.na(feature_matrix$ko_count)),
  stringsAsFactors = FALSE
)

cat("\n  Summary:\n")
for (i in 1:ncol(summary_stats)) {
  cat("    ", colnames(summary_stats)[i], ":", summary_stats[[i]], "\n")
}

# Save summary
fwrite(summary_stats, file.path(RESULTS_DIR, "06_summary_statistics.tsv"), sep = "\t")
cat("  ✓ Saved:", file.path(RESULTS_DIR, "06_summary_statistics.tsv"), "\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("\nAll outputs saved to:", RESULTS_DIR, "\n")
cat("\nNext steps:\n")
cat("  1. Install KEGGREST: BiocManager::install('KEGGREST')\n")
cat("  2. Run KofamScan on all protein sequences\n")
cat("  3. Map KOs to modules using KEGG Mapper or KEMET\n")
cat("  4. Add environment nutrient scores and carbon substrate diversity\n")
cat("  5. Proceed to modeling phase\n")
cat("================================================================================\n")

