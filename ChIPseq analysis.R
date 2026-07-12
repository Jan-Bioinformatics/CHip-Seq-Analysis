###############################################################################
# Differential ChIP-seq and Chromatin Accessibility Analysis at the MYC Locus: 
# Treatment-Induced TF Redistribution and Regulatory Remodelling

#                    Author: Janmita, Kaverimane Umesh
#                    Genome build: hg19

# Pipeline: QuasR/Rbowtie alignment → MACS3 peak calling → 
#           monaLisa/JASPAR2022 motif enrichment → ChIPseeker annotation 
#           → GREAT functional enrichment

###############################################################################

# Set directory
# setwd("C:/Users/Admin/OneDrive/Desktop/R studio Files/Applied Genomics/CHip Seq")

# 2. Load packages:
library(QuasR)          # Manages alignments and quality control tracking
library(Rsamtools)      # Processes low-level BAM/SAM sequencing files
library(GenomicAlignments)# Tracks and manipulates aligned genomic reads
library(rtracklayer)    # Imports/exports genomic data formats (BED, BigWig)
library(GenomicRanges)  # Handles genomic intervals, peak coordinates, and overlaps
library(BSgenome.Hsapiens.UCSC.hg19) # Accesses full hg19 human reference genome sequence
library(MACSr)          # Runs MACS3 peak-calling directly within R
library(dbplyr)         # Allows dplyr functions to work on database backends
library(BiocFileCache)  # Manages local caching for downloaded remote files
library(GenomeInfoDb)   # Standardises chromosome names across genome versions
library(monaLisa)       # Performs motif enrichment analysis and visualization
library(JASPAR2022)     # Comprehensive database of transcription factor profiles
library(TFBSTools)      # Analyses and manipulates transcription factor motifs
library(ChIPseeker)     # Annotates peaks relative to nearest genomic features
library(GenomicFeatures)# Creates and manipulates custom transcript databases
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # hg19 gene structure metadata (exons, promoters)
library(AnnotationDbi)  # Translates gene IDs (e.g., Entrez IDs to symbols)
library(ggrepel)        # Prevents text overlap in volcano plots and diagrams
library(ggplot2)        # Generates clean, publication-ready graphics and charts
library(org.Hs.eg.db)   # Maps Entrez IDs to symbols for human genes
library(patchwork)     # Combines multiple plots into a unified layout


################################################################################


################################################################################
##################### Step 1. OPTIONAL Fastq mapping. ##########################

# 1. Prepare FASTQ Sample File
#QuasR requires a **sample table**.
#samples <- data.frame( FileName = "test_ChIPseq.fastq.gz", SampleName = "ChIP_sample" )
#write.table(samples, file="samples.txt", sep="\t", row.names=FALSE, quote=FALSE)

#  2. Align FASTQ to hg19
#proj <- qAlign(sampleFile = "samples.txt", genome = "BSgenome.Hsapiens.UCSC.hg19", aligner = "Rbowtie", paired = "no")

# Output files: ChIP_sample.bam


#################################################################################
##################### Step 2.  Generating BigWig files ##############################
                  
# 1. Setup Directories 
base_path <- "C:/atac seq"
setwd(base_path)
dir.create("01_filtered_bam1", showWarnings = FALSE)
dir.create("02_bigwig1", showWarnings = FALSE)

samples <- c("ATACseq_treatment", "ATACseq_UNT")
genome_info <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)

for (sample in samples) {
  message("\n>>> Starting: ", sample)
  
  raw_input    <- paste0(sample, ".bam")
  filtered_bam <- file.path("01_filtered_bam", paste0(sample, "_filtered.bam"))
  dn_bam       <- file.path("01_filtered_bam", paste0(sample, "_downscaled.bam"))
  bw_output    <- file.path("02_bigwig", paste0(sample, "_normalized.bw"))
  
  # --- STEP 0: INDEXING ---
  if (!file.exists(paste0(raw_input, ".bai"))) {
    message("Indexing raw BAM...")
    indexBam(raw_input)
  }
  
  # --- STEP A: FILTERING (STREAMED) ---
  if (!file.exists(filtered_bam)) {
    message("Filtering MAPQ > 30...")
    filterBam(raw_input, destination = filtered_bam, 
              param = ScanBamParam(what = "mapq"),
              filter = FilterRules(list(high = function(x) x$mapq >= 30)))
  }
  
  # --- STEP B: DOWNSCALING ---
  if (!file.exists(dn_bam)) {
    message("Downscaling to 80%...")
    indexBam(filtered_bam)
    set.seed(123)
    filterBam(filtered_bam, destination = dn_bam, 
              filter = FilterRules(list(sub = function(x) runif(nrow(x)) < 0.80)))
    indexBam(dn_bam)
  }
  
  # --- STEP C: BIGWIG GENERATION (THE CRASH FIX) ---
  message("Generating BigWig (Memory-Safe Mode)...")
  
  # Instead of loading everything, we read it in and immediately convert to GRanges
  # We use 'use.names=FALSE' to save memory
  gr <- granges(readGAlignments(dn_bam, use.names=FALSE))
  
  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "tidy")
  
  # Match genome info
  seqlengths(gr) <- genome_info[seqlevels(gr)]
  
  # Calculate coverage and export IMMEDIATELY
  # We use a 'SimpleRleList' which is the most memory-efficient format
  cov_data <- coverage(gr)
  
  # Remove the GRanges object BEFORE exporting to free up those '48 bytes' and more
  rm(gr)
  gc()
  
  message("Exporting BigWig file...")
  export.bw(cov_data, con = bw_output)
  
  # --- STEP D: CLEANUP ---
  rm(cov_data)
  gc()
  message("Finished: ", sample)
}


#################################################################################
##################### Step 3. MAC peak calling ##############################


for (d in c("03_macs3_peaks", "04_merged_peaks", "05_counts",
            "06_differential", "07_motifs",
            "08_annotation", "09_plots")) {
  dir.create(d, showWarnings = FALSE)
}

# ============================================================
# 3. LOAD NARROWPEAK FILES
# ============================================================
# Note: Run on Windows; imported pre-computed MACS3 peak files to bypass macOS native requirements.

# Load the e6 peaks for each of the 4 TF ChIP samples
p_UNT1     <- import("03_macs3_peaks/UNT1_e8_peaks.narrowPeak")
p_UNT2     <- import("03_macs3_peaks/UNT2_e8_peaks.narrowPeak")
p_Treated1 <- import("03_macs3_peaks/treated_1_e8_peaks.narrowPeak")
p_Treated2 <- import("03_macs3_peaks/treated_2_e8_peaks.narrowPeak")

# How many peaks were called in each sample?
cat("=== Peaks per sample (e6 threshold) ===\n")
cat("UNT1:    ", length(p_UNT1),     "\n")
cat("UNT2:    ", length(p_UNT2),     "\n")
cat("Treated1:", length(p_Treated1), "\n")
cat("Treated2:", length(p_Treated2), "\n")

# ============================================================
# 4. MERGE PEAKS INTO ONE UNION SET 
# ============================================================
# Why: we need one consistent set of genomic regions so that
# every sample is compared on the same coordinates.
# reduce() merges any overlapping peaks into a single region.

all_peaks   <- c(p_UNT1, p_UNT2, p_Treated1, p_Treated2)
union_peaks <- reduce(all_peaks)

cat("\nTotal union peaks:", length(union_peaks), "\n")

export(union_peaks, "04_merged_peaks/union_peaks.bed")


# ============================================================
# 5. COUNT READS IN EACH PEAK FOR EACH SAMPLE  
# ============================================================
# We use the downscaled BAMs so library size is equal.
# countOverlaps() asks: for each peak, how many reads
# from this BAM file overlap with it?

genome_hg19 <- BSgenome.Hsapiens.UCSC.hg19

peaks <- import("04_merged_peaks/union_peaks.bed")

count_reads <- function(bam_path, peaks) {
  bam <- readGAlignments(bam_path)
  countOverlaps(peaks, bam)
}

cat("Counting reads — this may take a few minutes...\n")
counts_UNT1     <- count_reads("01_filtered_bam/UNT1_downscaled.bam",     peaks)
counts_UNT2     <- count_reads("01_filtered_bam/UNT2_downscaled.bam",     peaks)
counts_Treated1 <- count_reads("01_filtered_bam/Treated1_downscaled.bam", peaks)
counts_Treated2 <- count_reads("01_filtered_bam/Treated2_downscaled.bam", peaks)

# Build count table
count_df <- as.data.frame(peaks)[, 1:3]
colnames(count_df) <- c("chr", "start", "end")
count_df$UNT1     <- counts_UNT1
count_df$UNT2     <- counts_UNT2
count_df$Treated1 <- counts_Treated1
count_df$Treated2 <- counts_Treated2

write.table(count_df, "05_counts/peak_counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Q2a TABLE: peaks per sample and per group ----
summary_table1 <- data.frame(
  Sample = c("UNT1", "UNT2", "Treated1", "Treated2", 
             "UNT group (both)", "Treated group (both)"),
  
  N_peaks = c(
    sum(counts_UNT1 > 0),
    sum(counts_UNT2 > 0),
    sum(counts_Treated1 > 0),
    sum(counts_Treated2 > 0),
    # This counts rows where BOTH UNT1 and UNT2 are greater than 0
    sum(rowSums(count_df[, c("UNT1", "UNT2")] > 0) == 2),
    # This counts rows where BOTH Treated1 and Treated2 are greater than 0
    sum(rowSums(count_df[, c("Treated1", "Treated2")] > 0) == 2)
  )
)

print(summary_table1)

# Save the result
write.csv(summary_table1, "05_counts/peak_summary_table.csv", row.names = FALSE)



#################################################################################
##################### Step 4. DIFFERENTIAL PEAK ANALYSIS #######################


data <- read.table("05_counts/peak_counts.txt",
                   header = TRUE, sep = "\t")

data$mean_UNT     <- rowMeans(data[, c("UNT1","UNT2")])
data$mean_Treated <- rowMeans(data[, c("Treated1","Treated2")])

# Pseudocount of 1 avoids division by zero
data$fold_change <- (data$mean_Treated + 1) / (data$mean_UNT + 1)
data$log2FC      <- log2(data$fold_change)

# ---- FIXED t-test ----
# tryCatch wraps every t-test call.
# If the data are constant (identical values), it returns
# p = 1 instead of crashing. This is correct behaviour:
# if there is zero variance between replicates, there is
# no evidence of differential binding (p = 1).

data$pvalue <- apply(data, 1, function(x) {
  unt <- as.numeric(x[c("UNT1","UNT2")])
  trt <- as.numeric(x[c("Treated1","Treated2")])
  
  tryCatch(
    t.test(unt, trt, alternative = "two.sided")$p.value,
    error   = function(e) 1,   # constant data → p = 1
    warning = function(w) 1    # near-constant  → p = 1
  )
})

data$padj <- p.adjust(data$pvalue, method = "BH")

# Apply thresholds
diff_peaks <- subset(data,
                     (fold_change >= 2 | fold_change <= 0.5) & pvalue <= 0.05)

diff_peaks$direction <- ifelse(diff_peaks$log2FC > 0,
                               "gained", "lost")

write.table(diff_peaks, "06_differential/differential_peaks.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

n_total  <- nrow(diff_peaks)
n_gained <- sum(diff_peaks$direction == "gained")
n_lost   <- sum(diff_peaks$direction == "lost")

cat("\n=== Differential binding sites ===\n")
cat("Total:  ", n_total,  "\n")
cat("Gained: ", n_gained, "\n")
cat("Lost:   ", n_lost,   "\n")


# ============================================================
# 7. BARPLOT — GAINED VS LOST  
# ============================================================

plot_data <- data.frame(
  Direction = c("Gained", "Lost"),
  Count     = c(n_gained, n_lost)
)

p_bar <- ggplot(plot_data, aes(x = Direction, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.8, ) +
  geom_text(aes(label = Count), vjust = -0.4, size = 8) +
  scale_fill_manual(values = c("Gained" = "#1565C0",
                               "Lost"   = "#B71C1C")) +
  labs(
    title    = "Differential TF binding sites: treated vs untreated",
      x = NULL, y = "Number of peaks"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

ggsave("09_plots/Q2c_differential_barplot.pdf",
       p_bar, width = 3, height = 5)
print(p_bar)


#################################################

# 1. Start with your existing plot_data data frame
# (n_gained, n_lost, and n_total must already be defined in your environment)
# For this example, let's assume they are from your 129 vs 7 result.

plot_data <- data.frame(
  Direction = c("Gained", "Lost"),
  Count     = c(n_gained, n_lost)
)

# 2. Add Percentage and position data for labels (this makes it pro)
# 'ymax' and 'ymin' are required for ggplot's polar coordinate system
plot_data_pie <- plot_data %>%
  arrange(desc(Direction)) %>% # Standard order
  mutate(
    Percentage = (Count / sum(Count)) * 100,
    # Standard math for placing labels in the center of the slices
    LabelPos = cumsum(Count) - (0.5 * Count)
  )

# 3. Create the Publication-Quality Pie Chart
p_pie <- ggplot(plot_data_pie, aes(x = 1, y = Count, fill = Direction)) +
  # Make the base stacked bar
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
  # Transform to circular coordinates
  coord_polar(theta = "y") +
  
  # --- Add Labels: Both the Count and the Percentage ---
  # We use 'paste0' to combine the count with the (X.X%) text.
  # "Gained" will have "129 (94.9%)" and "Lost" will have "7 (5.1%)".
  geom_text(aes(y = LabelPos, 
                label = paste0(Count, "\n(", round(Percentage, 1), "%)")),
            size = 5, fontface = "bold", color = "white") +
  
  # --- Standard Publication Colors ---
  # Red for Gained (Increased Activity), Gray/Light Blue for Lost (Decreased Activity)
  scale_fill_manual(values = c("Gained" = "#e41a1c", # Distinct Red
                               "Lost"   = "#377eb8")) + # Distinct Blue
  
  # Clean, minimal theme suited for nature-styleComposition plots
  theme_void() + # This removes ALL lines, axes, and grids
  
  labs(
    title    = "Differential TF Binding Site Composition Treated Vs. Untreated",
    subtitle = paste0("Gain:FC ≥ 2 | Loss:FC ≤ 0.5  |  p ≤ 0.05  |  Total = ", n_total)
  ) +
  
  # Professional theme adjustments
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray30"),
    legend.position = "right", # Pie charts often need a legend for composition
    legend.title    = element_blank(),
    legend.text     = element_text(size = 15, color = "black")
  )

# 4. Save and Show
# PDFs are preferred for high-quality printing (vector format)
ggsave("09_plots/Q2c_differential_piechart.pdf", p_pie, width = 6, height = 5)
print(p_pie)



#################################################################################
##################### Step 5. MOTIF ANALYSIS  #######################


# Read narrowPeak
peaks_df <- read.table("06_differential/differential_peaks.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
peaks_df <- peaks_df[peaks_df$direction == "gained", ]
# MACS3 narrowPeak standard columns: 1=chr, 2=start, 3=end
colnames(peaks_df)[1:3] <- c("chrom", "start", "end")

# Create the GRanges object
gr_df <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns = TRUE, ignore.strand = TRUE)

# CRITICAL: Fix Chromosome Names (Changes '1' to 'chr1')
# This solves the "seqlengths" error by matching UCSC style
seqlevelsStyle(gr_df) <- "UCSC"

# Keep only standard chromosomes (1-22, X, Y, M) to remove random scaffolds
gr_df <- keepStandardChromosomes(gr_df, pruning.mode = "tidy")

# Prepare Sequences for Motif Search ---
# We center the search on 250bp around the peak middle (where TFs usually bind)
gr_centered <- resize(gr_df, width = 500, fix = "center")

#Sometimes extending peaks near chromosome ends creates invalid ranges.
#gr_centered <- trim(gr_centered)

# Extract the DNA sequences from the hg19 genome
# Explicitly set genome info to prevent sampling errors
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(gr_centered) <- seqlengths(genome)[seqlevels(gr_centered)]
peak_seqs <- getSeq(genome, gr_centered)

# confirm peak size
table(width(peak_seqs))

keep <- grepl("^[ACGT]+$", as.character(peak_seqs))
peak_seqs <- peak_seqs[keep]

# Since we want to analyze ALL peaks as one set, we create a factor 
# of the same length as our sequences, assigning them all to one bin.
my_bins <- factor(rep("all_peaks", length(peak_seqs)))

# Load Motif Database (JASPAR) ---
# We use the 'Vertebrates' collection from JASPAR
motifs <- getMatrixSet(JASPAR2022, opts = list(tax_group = "vertebrates", collection = "CORE"))
pwmL <- TFBSTools::toPWM(motifs)

# Run Enrichment ---
# Note the 'bins = my_bins' and 'background = "genome"'
# By choosing 'genome', monaLisa will compare your peaks 
# against random GC-matched regions in the genome.
motif_results <- calcBinnedMotifEnrR(
  seqs = peak_seqs,
  bins = my_bins,
  pwmL = pwmL,
  background = "genome",    # Best when using only one bin, # Changes 'genome' to 'shuffled'
  genome = BSgenome.Hsapiens.UCSC.hg19,
  #BPPARAM = SerialParam()
)

# Extract motif enrichment results:
log2_enrichment <- assay(motif_results, "log2enr")
head(log2_enrichment)

FDR <- assay(motif_results, "negLog10Padj")
head(FDR)

motif_info <- rowData(motif_results)
head(motif_info)

results_table_g <- data.frame(
  motif = rownames(motif_results),
  TF = motif_info$motif.name,
  log2_enrichment = apply(log2_enrichment, 1, max),
  minlog10_FDR = apply(FDR, 1, max)
)

head(results_table_g)


# Sort for top motifs
results_table_g <- results_table_g[order(-results_table_g$minlog10_FDR), ]
head(results_table_g, 20)

write.csv(results_table_g,"07_motifs/motifs_gained_genome.csv", row.names = FALSE)


# ============================================================
# LOST SITES ANALYSIS
# ============================================================

# filter for "lost" instead of "gained"
peaks_df_lost <- read.table("06_differential/differential_peaks.txt", 
                            sep="\t", header=TRUE, stringsAsFactors=FALSE)
peaks_df_lost <- peaks_df_lost[peaks_df_lost$direction == "lost", ]  

colnames(peaks_df_lost)[1:3] <- c("chrom", "start", "end")

# Check how many lost sites 
cat("Number of lost sites:", nrow(peaks_df_lost), "\n")
# Output: Number of lost sites: 7

#  Warning if too few peaks
if(nrow(peaks_df_lost) < 20){
  warning("Only ", nrow(peaks_df_lost), " lost sites found. 
           Motif results may not be statistically reliable.")
}

gr_df_lost <- makeGRangesFromDataFrame(peaks_df_lost,    
                                       keep.extra.columns=TRUE, 
                                       ignore.strand=TRUE)

seqlevelsStyle(gr_df_lost) <- "UCSC"                    
gr_df_lost <- keepStandardChromosomes(gr_df_lost, pruning.mode="tidy")

gr_centered_lost <- resize(gr_df_lost, width=500, fix="center")  

genome <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(gr_centered_lost) <- seqlengths(genome)[seqlevels(gr_centered_lost)]  

peak_seqs_lost <- getSeq(genome, gr_centered_lost)        

table(width(peak_seqs_lost))

keep_lost <- grepl("^[ACGT]+$", as.character(peak_seqs_lost))  
peak_seqs_lost <- peak_seqs_lost[keep_lost]

my_bins_lost <- factor(rep("all_peaks", length(peak_seqs_lost)))  


motif_results_lost <- calcBinnedMotifEnrR(               
  seqs = peak_seqs_lost,                                 
  bins = my_bins_lost,                                   
  pwmL = pwmL,
  background = "genome",
  genome = BSgenome.Hsapiens.UCSC.hg19
)

log2_enrichment_lost <- assay(motif_results_lost, "log2enr")   
FDR_lost <- assay(motif_results_lost, "negLog10Padj")          
motif_info_lost <- rowData(motif_results_lost)                 

results_table_lost <- data.frame( motif = rownames(motif_results_lost),
  TF = motif_info_lost$motif.name,
  log2_enrichment = apply(log2_enrichment_lost, 1, max),
  minlog10_FDR = apply(FDR_lost, 1, max)
)

results_table_lost <- results_table_lost[order(-results_table_lost$minlog10_FDR), ]
head(results_table_lost, 20)

# output filename
write.csv(results_table_lost, "07_motifs/motifs_lost_genome.csv", row.names=FALSE)

cat("\nNOTE: Lost site analysis based on only", nrow(peaks_df_lost), 
    "regions. Interpret results with caution.\n")



####################### Motif analysis when using total ATACseq peaks as background

# For gained peaks

p_ATAC_UNT     <- import("03_macs3_peaks/ATAC_UNT_e10_peaks.narrowPeak")
p_ATAC2_treated <- import("03_macs3_peaks/ATAC_treated_e10_peaks.narrowPeak")

library(GenomeInfoDb)

#Force consistent chromosome sets
p_ATAC_UNT <- keepStandardChromosomes(p_ATAC_UNT, pruning.mode="coarse")
p_ATAC2_treated <- keepStandardChromosomes(p_ATAC2_treated, pruning.mode="coarse")

#Ensure same naming style
seqlevelsStyle(p_ATAC_UNT) <- "UCSC"
seqlevelsStyle(p_ATAC2_treated) <- "UCSC"

#Merge the ATAC files
all_peaks <- reduce(c(p_ATAC_UNT, p_ATAC2_treated))

cat("\nTotal union peaks:", length(all_peaks), "\n")

export(all_peaks, "04_merged_peaks/all_peaks.bed")


# Step 1 — Load background peaks
bg <- read.table("04_merged_peaks/all_peaks.bed", sep="\t", header=TRUE)
colnames(bg)[1:3] <- c("chrom","start","end")

gr_bg <- makeGRangesFromDataFrame(bg, keep.extra.columns=TRUE, ignore.strand=TRUE)
seqlevelsStyle(gr_bg) <- "UCSC"
gr_bg <- keepStandardChromosomes(gr_bg, pruning.mode="tidy")

gr_bg_centered <- resize(gr_bg, width=500, fix="center")
gr_bg_centered <- trim(gr_bg_centered)                          # ← trim before getSeq

bg_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_bg_centered)
names(bg_seqs) <- as.character(gr_bg_centered)                  # ← assign names explicitly

keep_bg <- grepl("^[ACGT]+$", as.character(bg_seqs))
bg_seqs       <- bg_seqs[keep_bg]                               # ← keep DNAStringSet subset
gr_bg_centered <- gr_bg_centered[keep_bg]                       # ← keep GRanges in sync


# Step 2 — Load foreground (gained) peaks
peaks_df <- read.table("06_differential/differential_peaks.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
peaks_df <- peaks_df[peaks_df$direction == "gained", ]
colnames(peaks_df)[1:3] <- c("chrom", "start", "end")

gr_df <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns=TRUE, ignore.strand=TRUE)
seqlevelsStyle(gr_df) <- "UCSC"
gr_df <- keepStandardChromosomes(gr_df, pruning.mode="tidy")

gr_centered <- resize(gr_df, width=500, fix="center")
gr_centered <- trim(gr_centered)                                # ← trim before getSeq

peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_centered)
names(peak_seqs) <- as.character(gr_centered)                   # ← assign names explicitly

keep <- grepl("^[ACGT]+$", as.character(peak_seqs))
peak_seqs  <- peak_seqs[keep]                                   # ← keep DNAStringSet subset
gr_centered <- gr_centered[keep]                                # ← keep GRanges in sync


# Step 3 — Combine foreground + background
all_seqs <- c(peak_seqs, bg_seqs)

comb_bins <- factor(c(
  rep("foreground", length(peak_seqs)),
  rep("background", length(bg_seqs))
))

# Verify no empty names before passing to monaLisa
stopifnot(all(nchar(names(all_seqs)) > 0))
cat("Foreground:", length(peak_seqs), "| Background:", length(bg_seqs), "\n")

# step 4 - extract
motifs <- getMatrixSet(JASPAR2022, opts = list(tax_group = "vertebrates", collection = "CORE"))
pwmL <- TFBSTools::toPWM(motifs)


# Step 5 — Run monaLisa
motif_results <- calcBinnedMotifEnrR(
  seqs = all_seqs,
  bins = comb_bins,
  pwmL = pwmL
)

# Step 6 — Extract motif enrichment results:
log2_enrichment <- assay(motif_results, "log2enr")
head(log2_enrichment)

FDR <- assay(motif_results, "negLog10Padj")
head(FDR)

motif_info <- rowData(motif_results)
head(motif_info)

results_table <- data.frame(
  motif = rownames(motif_results),
  TF = motif_info$motif.name,
  log2_enrichment = apply(log2_enrichment, 1, max),
  minlog10_FDR = apply(FDR, 1, max)
)

head(results_table)


# Sort for top motifs
results_table <- results_table[order(-results_table$minlog10_FDR), ]
head(results_table, 20)

write.csv(results_table,"07_motifs/motifs_gained_ATAC.csv", row.names = FALSE)


# ============================================================
# LOST SITES — ATACseq Background
# Background (all_peaks) is already made — no need to redo!
# ============================================================

# Step 2 — Load foreground (LOST peaks)  
peaks_df_lost <- read.table("06_differential/differential_peaks.txt", 
                            sep="\t", header=TRUE, stringsAsFactors=FALSE)
peaks_df_lost <- peaks_df_lost[peaks_df_lost$direction == "lost", ]  
colnames(peaks_df_lost)[1:3] <- c("chrom", "start", "end")

gr_df_lost <- makeGRangesFromDataFrame(peaks_df_lost, keep.extra.columns=TRUE, ignore.strand=TRUE)
seqlevelsStyle(gr_df_lost) <- "UCSC"
gr_df_lost <- keepStandardChromosomes(gr_df_lost, pruning.mode="tidy")
gr_centered_lost <- resize(gr_df_lost, width=500, fix="center")
gr_centered_lost <- trim(gr_centered_lost)
peak_seqs_lost <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_centered_lost)
names(peak_seqs_lost) <- as.character(gr_centered_lost)
keep_lost <- grepl("^[ACGT]+$", as.character(peak_seqs_lost))
peak_seqs_lost   <- peak_seqs_lost[keep_lost]
gr_centered_lost <- gr_centered_lost[keep_lost]

cat("Lost foreground:", length(peak_seqs_lost), "| Background:", length(bg_seqs), "\n")

# Step 3 — Combine lost foreground + background
all_seqs_lost <- c(peak_seqs_lost, bg_seqs)                     
comb_bins_lost <- factor(c(
  rep("foreground", length(peak_seqs_lost)),                     
  rep("background", length(bg_seqs))
))

stopifnot(all(nchar(names(all_seqs_lost)) > 0))


# Step 5 — Run monaLisa
motif_results_lost <- calcBinnedMotifEnrR(                       
  seqs = all_seqs_lost,                                          
  bins = comb_bins_lost,                                         
  pwmL = pwmL
)

# Step 6 — Extract results
log2_enrichment_lost <- assay(motif_results_lost, "log2enr")
FDR_lost             <- assay(motif_results_lost, "negLog10Padj")
motif_info_lost      <- rowData(motif_results_lost)

results_table_lost <- data.frame(
  motif           = rownames(motif_results_lost),
  TF              = motif_info_lost$motif.name,
  log2_enrichment = apply(log2_enrichment_lost, 1, max),
  minlog10_FDR    = apply(FDR_lost, 1, max)
)

# Sort by significance
results_table_lost <- results_table_lost[order(-results_table_lost$minlog10_FDR), ]
head(results_table_lost, 20)

write.csv(results_table_lost, "07_motifs/motifs_lost_ATAC.csv", row.names=FALSE)  


# ============================================================
#       VOLCANO PLOT: -log10(FDR) vs log(Enrichment)
# ============================================================

# --- Plot 1: Genome background results ---
results_table_g$significant <- results_table_g$minlog10_FDR > 2
# Convert log2 → natural log by multiplying by ln(2)
results_table_g$log_enrichment <- results_table_g$log2_enrichment * log(2)
results_table$log_enrichment   <- results_table$log2_enrichment   * log(2)
top_labels_g <- head(results_table_g[order(-results_table_g$minlog10_FDR), ], 10)

p1 <- ggplot(results_table_g, aes(x = log_enrichment, y = minlog10_FDR)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#D85A30", "FALSE" = "#B4B2A9"),
                     labels = c("TRUE" = "FDR < 0.01", "FALSE" = "Not significant")) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "#888780", linewidth = 0.5) +
  geom_text_repel(data = top_labels_g,
                  aes(label = TF),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.4,
                  segment.color = "#888780",
                  segment.size = 0.3) +
  labs(title = "Motif Enrichment – Genome Background",
       subtitle = "Gained ATAC-seq peaks vs random genomic regions",
       x = "log(Enrichment)",                         
       y = "-log10(adjusted p-value)",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

print(p1)
ggsave("07_motifs/volcano_gained_genome.pdf", p1, width = 7, height = 6)

# --- Plot 2: ATAC background results ---
results_table$significant <- results_table$minlog10_FDR > 2
results_table_g$log_enrichment <- results_table_g$log2_enrichment * log(2)     # natural log
top_labels_a <- head(results_table[order(-results_table$minlog10_FDR), ], 10)

p2 <- ggplot(results_table, aes(x = log_enrichment, y = minlog10_FDR)) +
  geom_point(aes(color = significant), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("TRUE" = "#185FA5", "FALSE" = "#B4B2A9"),
                     labels = c("TRUE" = "FDR < 0.01", "FALSE" = "Not significant")) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "#888780", linewidth = 0.5) +
  geom_text_repel(data = top_labels_a,
                  aes(label = TF),
                  size = 3,
                  max.overlaps = 15,
                  box.padding = 0.4,
                  segment.color = "#888780",
                  segment.size = 0.3) +
  labs(title = "Motif Enrichment – ATAC Background",
       subtitle = "Gained ATAC-seq peaks vs total ATAC-seq union peaks",
       x = "log(Enrichment)",                        
       y = "-log10(adjusted p-value)",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

print(p2)
ggsave("07_motifs/volcano_gained_ATAC.pdf", p2, width = 7, height = 6)

# --- Combined side-by-side panel ---
combined <- p1 + p2 +
  plot_annotation(title = "Motif Enrichment at Gained Binding Sites",
                  tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("07_motifs/volcano_combined.pdf", combined, width = 14, height = 6)



#################################################################################
##################### Step 6. ANNOTATE DIFFERENTIAL BINDING SITES  #######################


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# ============================================================
# STEP 1: Prepare GRanges for gained and lost peaks separately
# ============================================================

gained_peaks <- diff_peaks[diff_peaks$direction == "gained", ]
lost_peaks   <- diff_peaks[diff_peaks$direction == "lost",   ]

gr_gained <- makeGRangesFromDataFrame(gained_peaks,
                                      seqnames.field = "chr",
                                      start.field    = "start",
                                      end.field      = "end",
                                      keep.extra.columns = TRUE)

gr_lost <- makeGRangesFromDataFrame(lost_peaks,
                                    seqnames.field = "chr",
                                    start.field    = "start",
                                    end.field      = "end",
                                    keep.extra.columns = TRUE)

seqlevelsStyle(gr_gained) <- "UCSC"
seqlevelsStyle(gr_lost)   <- "UCSC"

gr_gained <- keepStandardChromosomes(gr_gained, pruning.mode = "tidy")
gr_lost   <- keepStandardChromosomes(gr_lost,   pruning.mode = "tidy")

# ============================================================
# STEP 2: Annotate with ChIPseeker
# ============================================================

anno_gained <- annotatePeak(gr_gained,
                            tssRegion  = c(-3000, 3000),
                            TxDb       = txdb,
                            annoDb     = "org.Hs.eg.db")

anno_lost   <- annotatePeak(gr_lost,
                            tssRegion  = c(-3000, 3000),
                            TxDb       = txdb,
                            annoDb     = "org.Hs.eg.db")

# Save annotated tables
write.csv(as.data.frame(anno_gained), "08_annotation/annotation_gained.csv", row.names = FALSE)
write.csv(as.data.frame(anno_lost),   "08_annotation/annotation_lost.csv",   row.names = FALSE)

# ============================================================
# STEP 3: Extract and clean annotation categories
# ============================================================

clean_annotation <- function(anno_df) {
  anno_df$annotation_clean <- gsub(" \\(.*\\)", "", anno_df$annotation)
  anno_df$annotation_clean <- dplyr::case_when(
    grepl("Promoter",   anno_df$annotation_clean) ~ "Promoter",
    grepl("Exon",       anno_df$annotation_clean) ~ "Exon",
    grepl("Intron",     anno_df$annotation_clean) ~ "Intron",
    grepl("Downstream", anno_df$annotation_clean) ~ "Downstream",
    grepl("Distal",     anno_df$annotation_clean) ~ "Distal Intergenic",
    grepl("3' UTR",     anno_df$annotation_clean) ~ "3' UTR",
    grepl("5' UTR",     anno_df$annotation_clean) ~ "5' UTR",
    TRUE ~ "Distal Intergenic"
  )
  return(anno_df)
}

df_gained <- clean_annotation(as.data.frame(anno_gained))
df_lost   <- clean_annotation(as.data.frame(anno_lost))

# ============================================================
# STEP 4: Calculate proportions for plotting
# ============================================================

get_proportions <- function(df, group_label) {
  tbl <- as.data.frame(table(df$annotation_clean))
  colnames(tbl) <- c("Feature", "Count")
  tbl$Percentage <- tbl$Count / sum(tbl$Count) * 100
  tbl$Group <- group_label
  tbl <- tbl[order(-tbl$Count), ]
  return(tbl)
}

prop_gained <- get_proportions(df_gained, "Gained")
prop_lost   <- get_proportions(df_lost,   "Lost")
prop_all    <- rbind(prop_gained, prop_lost)

# ============================================================
# STEP 5: BARPLOT — side by side
# ============================================================

# Consistent colour palette across all plots
feature_colors <- c(
  "Promoter"          = "#E63946",
  "Intron"            = "#457B9D",
  "Exon"              = "#2A9D8F",
  "Distal Intergenic" = "#E9C46A",
  "Downstream"        = "#F4A261",
  "3' UTR"            = "#A8DADC",
  "5' UTR"            = "#6A994E"
)

p_bar_anno <- ggplot(prop_all, aes(x = reorder(Feature, -Percentage),
                                   y = Percentage,
                                   fill = Feature)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            vjust = -0.4, size = 3.2) +
  facet_wrap(~Group, ncol = 2) +
  scale_fill_manual(values = feature_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Genomic Distribution of Differential Binding Sites",
       subtitle = "Gained vs Lost ATAC-seq peaks (treated vs untreated)",
       x = NULL,
       y = "Percentage of peaks (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 35, hjust = 1),
        legend.position = "none")

print(p_bar_anno)
ggsave("08_annotation/Q4a_barplot_annotation.pdf",
       p_bar_anno, width = 10, height = 6)

# ============================================================
# STEP 6: PIE CHARTS — gained and lost side by side
# ============================================================

make_pie <- function(prop_df, title_str) {
  prop_df <- prop_df[order(-prop_df$Count), ]
  prop_df$label_pos <- cumsum(prop_df$Percentage) - 0.5 * prop_df$Percentage
  
  ggplot(prop_df, aes(x = 1, y = Percentage, fill = Feature)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.4) +
    coord_polar(theta = "y") +
    geom_text(aes(y = label_pos,
                  label = ifelse(Percentage > 4,
                                 paste0(round(Percentage, 1), "%"), "")),
              size = 3.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = feature_colors) +
    theme_void(base_size = 11) +
    labs(title = title_str, fill = "Genomic Feature") +
    theme(plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.title = element_text(face = "bold"))
}

pie_gained <- make_pie(prop_gained, "Gained Peaks")
pie_lost   <- make_pie(prop_lost,   "Lost Peaks")

combined_pie <- pie_gained + pie_lost +
  plot_annotation(
    title    = "Genomic Feature Distribution of Differential Binding Sites",
    subtitle = "Gained vs Lost ATAC-seq peaks (treated vs untreated)",
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(face = "bold"))

print(combined_pie)
ggsave("08_annotation/Q4a_piechart_annotation.pdf",
       combined_pie, width = 12, height = 6)

# ============================================================
# STEP 7: SUMMARY TABLE — print to console
# ============================================================

cat("\n=== Genomic Distribution Summary ===\n\n")
cat("--- GAINED peaks ---\n")
print(prop_gained[, c("Feature", "Count", "Percentage")], row.names = FALSE)

cat("\n--- LOST peaks ---\n")
print(prop_lost[, c("Feature", "Count", "Percentage")], row.names = FALSE)


# ============================================================
#  ANNOTATE DIFFERENTIAL BINDING SITES (ALL PEAKS COMBINED)
# ============================================================

# STEP 1: Prepare GRanges for ALL differential peaks

gr_diff <- makeGRangesFromDataFrame(diff_peaks,
                                    seqnames.field = "chr",
                                    start.field    = "start",
                                    end.field      = "end",
                                    keep.extra.columns = TRUE)

seqlevelsStyle(gr_diff) <- "UCSC"
gr_diff <- keepStandardChromosomes(gr_diff, pruning.mode = "tidy")


# STEP 2: Annotate with ChIPseeker

anno_diff <- annotatePeak(gr_diff,
                          tssRegion = c(-3000, 3000),
                          TxDb      = txdb,
                          annoDb    = "org.Hs.eg.db")

# Save annotated table
write.csv(as.data.frame(anno_diff),
          "08_annotation/combined_annotation_all_diff_peaks.csv",
          row.names = FALSE)


# STEP 3: Clean annotation categories

df_diff <- as.data.frame(anno_diff)

df_diff$annotation_clean <- dplyr::case_when(
  grepl("Promoter",   df_diff$annotation) ~ "Promoter",
  grepl("5' UTR",     df_diff$annotation) ~ "5' UTR",
  grepl("3' UTR",     df_diff$annotation) ~ "3' UTR",
  grepl("Exon",       df_diff$annotation) ~ "Exon",
  grepl("Intron",     df_diff$annotation) ~ "Intron",
  grepl("Downstream", df_diff$annotation) ~ "Downstream",
  TRUE                                    ~ "Distal Intergenic"
)


# STEP 4: Calculate proportions

prop_df <- as.data.frame(table(df_diff$annotation_clean))
colnames(prop_df) <- c("Feature", "Count")
prop_df$Percentage <- prop_df$Count / sum(prop_df$Count) * 100
prop_df <- prop_df[order(-prop_df$Count), ]

# Print summary
cat("\n=== Genomic Distribution of All Differential Peaks ===\n")
print(prop_df, row.names = FALSE)

# ============================================================
# STEP 5: BARPLOT
# ============================================================

feature_colors <- c(
  "Promoter"          = "#E63946",
  "Intron"            = "#457B9D",
  "Exon"              = "#2A9D8F",
  "Distal Intergenic" = "#E9C46A",
  "Downstream"        = "#F4A261",
  "3' UTR"            = "#A8DADC",
  "5' UTR"            = "#6A994E"
)

p_bar <- ggplot(prop_df, aes(x = reorder(Feature, -Percentage),
                             y = Percentage,
                             fill = Feature)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),
            vjust = -0.4, size = 4) +
  scale_fill_manual(values = feature_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Genomic Distribution of Differential Binding Sites",
       subtitle = paste0("All differential peaks (n = ", nrow(df_diff), ")"),
       x = NULL,
       y = "Percentage of peaks (%)") +
  theme_classic(base_size = 13) +
  theme(plot.title      = element_text(face = "bold"),
        axis.text.x     = element_text(angle = 35, hjust = 1),
        legend.position = "none")

print(p_bar)
ggsave("08_annotation/Q4_Combined_barplot_annotation.pdf",
       p_bar, width = 7, height = 6)

# ============================================================
# STEP 6: PIE CHART
# ============================================================

prop_df$label_pos <- cumsum(prop_df$Percentage) - 0.5 * prop_df$Percentage

p_pie <- ggplot(prop_df, aes(x = 1, y = Percentage, fill = Feature)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.4) +
  coord_polar(theta = "y") +
  geom_text(aes(y = label_pos,
                label = ifelse(Percentage > 4,
                               paste0(round(Percentage, 1), "%"), "")),
            size = 4, fontface = "bold", color = "white") +
  scale_fill_manual(values = feature_colors) +
  theme_void(base_size = 12) +
  labs(title    = "Genomic Distribution of Differential Binding Sites",
       subtitle = paste0("All differential peaks (n = ", nrow(df_diff), ")"),
       fill     = "Genomic Feature") +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
        legend.title  = element_text(face = "bold"))

print(p_pie)
ggsave("08_annotation/Q4Combined_piechart_annotation.pdf",
       p_pie, width = 7, height = 6)




#################################################################################
##################### Step 7. GREAT ANALYSIS ##################################

# Fix chromosome names to UCSC style (chr1, chr2...) — all differential peaks combined
great_bed <- diff_peaks[, c("chr", "start", "end")]

# Add "chr" prefix if not already present
great_bed$chr <- ifelse(grepl("^chr", great_bed$chr),
                        great_bed$chr,
                        paste0("chr", great_bed$chr))

write.table(great_bed,
            "08_annotation/combined_great_input_all_diff_peaks.bed",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Verify
cat("Total peaks exported:", nrow(great_bed), "\n")
head(great_bed)