library(ape)
library(pegas)
library(readr)

# ==========================================================
# --- 1. SETTINGS & FOLDER SELECTION ---
# ==========================================================
PROJECT_LABEL <- "ITS1_5.8S_ITS2"  
QC_THRESHOLD  <- 0.45    
IS_CODING     <- FALSE   # TRUE for COI, FALSE for 18S/ITS, For partition of ITS--> ITS1_5.8S_ITS2
GENETIC_CODE  <- 5       # 5 = Invertebrate Mito, 1 = Standard Nuclear

cat("SELECT THE FOLDER WHERE YOU WANT TO SAVE ALL YOUR FILES...")
my_folder <- choose.dir(caption = "Select Results Folder")
setwd(my_folder)

#To read the fasta file properly

load_clean_mat <- function(path) {
  lines <- readLines(path)
  h_idx <- grep("^>", lines)
  names <- gsub("^>", "", lines[h_idx])
  #Strip everything from the "/" to the end (Jalview length suffix)
  names <- gsub("/.*$", "", names)
  seqs <- sapply(1:length(h_idx), function(i) {
    start <- h_idx[i] + 1
    end <- if(i < length(h_idx)) h_idx[i+1] - 1 else length(lines)
    paste(lines[start:end], collapse = "")
  })
  seqs <- toupper(seqs)
  seqs <- gsub("[^ACGT]", "-", seqs)
  mat <- do.call(rbind, strsplit(seqs, ""))
  rownames(mat) <- names
  return(mat)
}

# ==========================================================
# --- 2. STEP 1: INITIAL QC ---
# ==========================================================
cat("\nStep 1: Select your RAW UNFILTERED alignment generated from the Python API script...")

#Fasta file 
raw_path <- file.choose()
mat_raw  <- load_clean_mat(raw_path)

cat("\nStep 2: Select your raw data generated from the Python API script...")
#Metadata file
metadata_raw <- read_tsv(file.choose())

#Check that they are in the same order
if (!all(rownames(mat_raw) == metadata_raw$Sequence_ID)) {
  stop("The DNA and Metadata are not in the same order!")
}

#For gap percentage
total_start <- nrow(mat_raw)
gap_perc <- rowSums(mat_raw == "-") / ncol(mat_raw)
is_fragment <- gap_perc > QC_THRESHOLD

# Deletion Log
deletion_log <- data.frame(
  Sequence_ID = rownames(mat_raw),
  Country = metadata_raw$Country,
  Status      = ifelse(is_fragment, "Deleted", "Kept"),
  Reason      = ifelse(is_fragment, "Fragment (>45% gaps)", "Passed QC"),
  Gap_Percent = round(gap_perc * 100, 2),
  stringsAsFactors = FALSE
)
write.table(deletion_log, paste0(PROJECT_LABEL, "_Deletion_Log.tsv"), sep="\t", row.names=F, quote=F)

mat_passed_qc <- mat_raw[!is_fragment, , drop = FALSE]
metadata_passed_qc <- metadata_raw[!is_fragment, , drop = FALSE]


#Check that the metadata and fasta have the same length
nrow(mat_passed_qc) 
nrow(metadata_passed_qc)


# ==========================================================
# --- 3. STEP 2: PRE-TRIM COLLAPSE ---
# ==========================================================
#This is the step where we loose all original names aka accession numbers and get the Hap_1 names

haps_pre <- haplotype(as.DNAbin(mat_passed_qc), strict = FALSE)
idx_pre  <- attr(haps_pre, "index") #the map of who's identical
num_unique_pre <- nrow(haps_pre)

#Gives the first of the list of collapsed haplotypes, so the fasta will only have the representative's names, x[1]
pre_trim_rep_rows <- sapply(idx_pre, function(x) x[1])
mat_for_trimal <- mat_passed_qc[pre_trim_rep_rows, , drop = FALSE]

#Keeps the metadata of the representative haplotypes only
metadata_pre <- metadata_passed_qc[pre_trim_rep_rows, , drop = FALSE]
write.table(metadata_pre, paste0(PROJECT_LABEL, "_metadata_for_trimal.tsv"), sep="\t", row.names=F, quote=F)

pre_trim_fasta <- paste0(PROJECT_LABEL, "_for_trimal.fasta")
writeLines(paste0(">", rownames(mat_for_trimal), "\n", apply(mat_for_trimal, 1, paste, collapse = "")), pre_trim_fasta)

cat("\n--- ACTION REQUIRED ---")
cat("\n1. Run trimAl on:", pre_trim_fasta)
cat("\n2. Save it, then come back here to R.")

# ==========================================================
# --- 4. STEP 3: POST-TRIM & AUTO-FRAME FIXER ---
# ==========================================================
readline(prompt="Press [Enter] AFTER saving your trimmed file to continue...")

cat("\nStep 3: Select your TRIMMED alignment...")
trim_path <- file.choose()
mat_trim  <- load_clean_mat(trim_path)

final_offset <- 0
final_remainder <- 0

if (IS_CODING) {
  cat("\n--- ANALYZING READING FRAME ---")
  count_stops <- function(mat, offset, code) {
    test_mat <- mat
    if(offset > 0) test_mat <- test_mat[, (offset + 1):ncol(test_mat), drop = FALSE]
    rem <- ncol(test_mat) %% 3
    if(rem > 0) test_mat <- test_mat[, 1:(ncol(test_mat) - rem), drop = FALSE]
    prot <- trans(as.DNAbin(test_mat), code = code)
    prot_char <- as.character(prot)
    sum(apply(prot_char, 1, function(x) any(x[1:(length(x)-1)] == "*")))
  }
  
  s1 <- count_stops(mat_trim, 0, GENETIC_CODE)
  s2 <- count_stops(mat_trim, 1, GENETIC_CODE)
  s3 <- count_stops(mat_trim, 2, GENETIC_CODE)
  cat("\nStops detected: F1:", s1, "| F2:", s2, "| F3:", s3)
  
  final_offset <- which.min(c(s1, s2, s3)) - 1
  if (final_offset > 0) mat_trim <- mat_trim[, (final_offset + 1):ncol(mat_trim), drop = FALSE]
  
  final_remainder <- ncol(mat_trim) %% 3
  if (final_remainder != 0) mat_trim <- mat_trim[, 1:(ncol(mat_trim) - final_remainder), drop = FALSE]
}

# Align names back to current representatives
rownames(mat_trim) <- rownames(mat_for_trimal)

# Final Haplotype Collapse
haps_post <- haplotype(as.DNAbin(mat_trim), strict = FALSE)
idx_post  <- attr(haps_post, "index") # A list of numbers


# ==========================================================
# --- 5. MAPPING OF REPRESENTATIVE HAPLOTYPES ---
# ==========================================================
receipt_list  <- list()
long_map_list <- list()
supp_map_list <- list()
cumulative_n  <- c()

for (i in 1:length(idx_post)) {
  merged_rep_indices <- idx_post[[i]]
  rep_name <- rownames(mat_for_trimal)[merged_rep_indices[1]]
  all_original_names <- c()
  for (idx in merged_rep_indices) {
    original_names <- rownames(mat_passed_qc)[idx_pre[[idx]]]
    all_original_names <- c(all_original_names, original_names)
  }
  
  # 1. Grab all countries for this haplotype's original names
  group_countries <- metadata_passed_qc$Country[metadata_passed_qc$Sequence_ID %in% all_original_names]
  
  # 2. Clean the list: remove NAs and remove duplicates
  unique_countries <- unique(na.omit(group_countries))
  
  # 3. Glue them into one clean string (e.g., "Greece, Italy")
  final_country_string <- paste(unique_countries, collapse = ", ")

  cumulative_n[i] <- length(all_original_names)
  receipt_list[[i]] <- data.frame(Final_Hap_ID = paste0("Hap_", i),
                                  Total_n = length(all_original_names),
                                  Countries = final_country_string,
                                  Original_Names = paste(all_original_names, collapse = "; "),
                                  stringsAsFactors = FALSE)
  long_map_list[[i]] <- data.frame(Original_Sequence = all_original_names,
                                   Final_Hap_ID = paste0("Hap_", i),
                                   stringsAsFactors = FALSE)


# Supplementary Map (Representative vs. Collapsed)
  collapsed_others <- all_original_names[all_original_names != rep_name]
  if(length(collapsed_others) == 0) {collapsed_others <- "None (Unique)"
  
  }
  
  supp_map_list[[i]] <- data.frame(
    Haplotype_ID   = paste0("Hap_", i),
    Representative = rep_name,
    Collapsed_IDs  = paste(collapsed_others, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

master_receipt <- do.call(rbind, receipt_list)
supp_map       <- do.call(rbind, supp_map_list)

write.table(master_receipt, paste0(PROJECT_LABEL, "_final_metadata.tsv"), sep="\t", row.names=F, quote=F)
write.table(supp_map, paste0(PROJECT_LABEL, "_representative_map.tsv"), sep="\t", row.names=F, quote=F)

# Final Fasta
unique_rows <- sapply(idx_post, function(x) x[1])
final_mat <- mat_trim[unique_rows, , drop = FALSE]
rownames(final_mat) <- paste0("Hap_", 1:nrow(final_mat), "_n", cumulative_n)
writeLines(paste0(">", rownames(final_mat), "\n", apply(final_mat, 1, paste, collapse = "")), 
           paste0(PROJECT_LABEL, "_final_tree.fasta"))


final_summary_table <- data.frame(
  Final_Hap_ID   = paste0("Hap_", 1:length(receipt_list), "_n", cumulative_n),
  acc_versioned  = sapply(supp_map_list, function(x) x$Representative),
  species_name   = sapply(supp_map_list, function(x) {
                     metadata_passed_qc$Species[metadata_passed_qc$Sequence_ID == x$Representative][1]
                   }),
  Total_n        = cumulative_n,
  ncbi_country   = sapply(receipt_list, function(x) x$Countries),
  Original_IDs   = sapply(receipt_list, function(x) x$Original_Names),
  stringsAsFactors = FALSE
)

# Save this as your "Plotting Metadata"
write_tsv(final_summary_table, paste0(PROJECT_LABEL, "_Metadata_for_Plotting.tsv"))

# ==========================================================
# --- 6. MERGE TRACKER ---
# ==========================================================
merge_tracker <- data.frame(Original_Sequence = rownames(mat_passed_qc), Step_1_Initial_Group = NA,
                            Step_2_Post_Trim_Group = NA, Final_Hap_ID = NA, stringsAsFactors = FALSE)

for (i in seq_along(idx_pre)) {
  member_names <- rownames(mat_passed_qc)[idx_pre[[i]]]
  merge_tracker$Step_1_Initial_Group[merge_tracker$Original_Sequence %in% member_names] <- paste0("PreTrim_Group_", i)
}

for (i in seq_along(idx_post)) {
  merged_reps <- rownames(mat_trim)[idx_post[[i]]]
  all_members <- c()
  for (r_name in merged_reps) {
    group_idx <- which(sapply(idx_pre, function(x) rownames(mat_passed_qc)[x[1]] == r_name))
    all_members <- c(all_members, rownames(mat_passed_qc)[idx_pre[[group_idx]]])
  }
  merge_tracker$Step_2_Post_Trim_Group[merge_tracker$Original_Sequence %in% all_members] <- paste0("PostTrim_Group_", i)
  merge_tracker$Final_Hap_ID[merge_tracker$Original_Sequence %in% all_members] <- paste0("Hap_", i)
}

merge_tracker$Merge_Type <- "Unique"
dup_pre <- merge_tracker$Step_1_Initial_Group[duplicated(merge_tracker$Step_1_Initial_Group)]
merge_tracker$Merge_Type[merge_tracker$Step_1_Initial_Group %in% dup_pre] <- "Initial_Identity"

for (post_grp in unique(merge_tracker$Step_2_Post_Trim_Group)) {
  unique_pre_in_post <- unique(merge_tracker$Step_1_Initial_Group[merge_tracker$Step_2_Post_Trim_Group == post_grp])
  if (length(unique_pre_in_post) > 1) merge_tracker$Merge_Type[merge_tracker$Step_2_Post_Trim_Group == post_grp] <- "Merged_by_Trimming"
}
write.table(merge_tracker, paste0(PROJECT_LABEL, "_Merge_Event_Log.tsv"), sep="\t", row.names=F, quote=F)

# ==========================================================
# --- 7. COMPREHENSIVE FINAL SUMMARY & AUTO-CLEAN ---
# ==========================================================
verification <- "PASSED"
has_stop <- rep(FALSE, nrow(final_mat)) # Default: no stops

if(IS_CODING) {
  final_prot <- trans(as.DNAbin(final_mat), code = GENETIC_CODE)
  # Check each sequence for internal stops
  has_stop <- apply(as.character(final_prot), 1, function(x) any(x[1:(length(x)-1)] == "*"))
  
  if(any(has_stop)) {
    verification <- "FAILED - STOPS REMOVED"
    bad_hap_ids <- rownames(final_mat)[has_stop]
    cat("\n🚨 CULPRITS FOUND! Removing haplotypes with stop codons:", paste(bad_hap_ids, collapse=", "))
    
    # Actually remove them
    final_mat <- final_mat[!has_stop, , drop = FALSE]
  }
} else {
  verification <- "N/A (Non-coding)"
}

summary_text <- capture.output({
    cat("\n--- FINAL PIPELINE REPORT ---\n")
    cat("Status (Reading Frame):  ", verification, "\n")
    cat("1. Total Started:         ", total_start, "\n")
    cat("2. Passed QC (>45% gaps):", nrow(mat_passed_qc), "\n")
    cat("3. Haplotypes Pre-Trim:   ", num_unique_pre, "\n")
    cat("4. Haplotypes Post-Trim: ", length(has_stop), "\n")
    cat("5. Sequences Removed:    ", sum(has_stop), "\n")
    cat("6. Final Tree Sequences: ", nrow(final_mat), "\n")
    
    if(IS_CODING) {
      cat("7. Left-Side Shift:       ", final_offset, "bp\n")
      cat("8. Right-Side Trim:       ", final_remainder, "bp\n")
    }
    cat("------------------------------\n")
})

# Print to console
cat(summary_text, sep="\n")

# SAVE THE CLEANED FASTA
writeLines(paste0(">", rownames(final_mat), "\n", apply(final_mat, 1, paste, collapse = "")), 
           paste0(PROJECT_LABEL, "_final_tree_CLEANED.fasta"))

# Save summary
writeLines(summary_text, paste0(PROJECT_LABEL, "_Final_Summary.txt"))

# ==========================================================
# --- 8. IQ-TREE PARTITION GENERATOR (ADAPTIVE) ---
# ==========================================================

final_len <- ncol(final_mat)
partition_filename <- paste0(PROJECT_LABEL, "_partition.nex")

if (IS_CODING) {
  # --- COI / PROTEIN CODING LOGIC ---
  if (final_len %% 3 == 0) {
    partition_content <- c(
      "#nexus",
      "begin sets;",
      paste0("    charset part1 = 1-", final_len, "\\3;"),
      paste0("    charset part2 = 2-", final_len, "\\3;"),
      paste0("    charset part3 = 3-", final_len, "\\3;"),
      "end;"
    )
    writeLines(partition_content, partition_filename)
    cat("\nSUCCESS: NEXUS Codon Partition created:", partition_filename)
  } else {
    cat("\nERROR: Length", final_len, "is not divisible by 3.")
  }
  
} else if (PROJECT_LABEL == "ITS1_5.8S_ITS2") {
  # --- ITS (SEGMENTED) LOGIC ---
  cat("\n--- ITS PARTITION SETUP ---")
  cat("\nYour total alignment length is:", final_len)
  
  # We use s58_end instead of 5.8S_end because R variables can't start with numbers
  val_its1 <- as.numeric(readline(prompt="Enter the LAST column of ITS1 (e.g. 250): "))
  val_s58  <- as.numeric(readline(prompt="Enter the LAST column of 5.8S (e.g. 410): "))
  
  partition_content <- c(
    "#nexus",
    "begin sets;",
    paste0("    charset ITS1 = 1-", val_its1, ";"),
    paste0("    charset S5.8 = ", val_its1 + 1, "-", val_s58, ";"),
    paste0("    charset ITS2 = ", val_s58 + 1, "-", final_len, ";"),
    "end;"
  )
  writeLines(partition_content, partition_filename)
  cat("\nSUCCESS: NEXUS ITS Segmented Partition created:", partition_filename)

} else {
  # --- NO NEXUS CREATED ---
  cat("\nNOTICE: No partition file required for", PROJECT_LABEL, ". Skipping Section 8.")
}