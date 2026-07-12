# End-to-End ChIP-seq Analysis Workflow: Transcription Factor Dynamics

This repository contains a professional bioinformatics pipeline analyzing Transcription Factor (TF) binding dynamics under treated vs. untreated conditions. Using R/Bioconductor and MACS3, the project identifies high-confidence binding landscapes, measures quantitative differential binding shifts, discovers driver motifs, and maps downstream functional pathways.

---

## 📊 Pipeline Roadmap
1. **Pre-processing**: Quality filtering (MAPQ ≥ 30) and sequencing depth normalization via random downsampling.
2. **Peak Calling**: Positional tracking of binding landscapes using MACS3.
3. **Differential Metrics**: Inter-group quantitative analysis identifying high-confidence gained and lost binding domains.
4. **Motif Discovery**: GC-matched genomic and ATAC-seq background comparisons using `monaLisa`.
5. **Functional Annotation**: Genomic context mapping and pathway enrichment via `ChIPseeker` and `GREAT`.

---

## 🔬 Core Findings & Visualizations

### 1. Unified Locus Tracking (The MYC Locus)
Active transcriptional regulation and chromatin accessibility profiling are benchmarked across the *MYC* proto-oncogene locus. The interaction between accessible chromatin (`ATAC-seq`) and partial repressive marks (`H3K27me3`) highlights an epigenetically poised transcriptional landscape.

![MYC Locus IGV Snapshot](09_plots/myc_locus_igv.png)
*Figure 1: Comprehensive IGV tracking across the MYC locus validating ChIP-seq, ATAC-seq, and ENCODE histone signals.*

### 2. Quantitative Inter-Group Variations
Using strict exploratory filters ($p \le 0.05$; $\text{Fold Change} \ge 2$), a comprehensive comparison reveals a massive quantitative shift toward open chromatin and active binding configurations post-treatment:
- **Total Selected Peaks**: 14,170
- **Statistically Significant Shifts**: 136 sites
- **Gained Sites (94.9%)**: 129 peaks
- **Lost Sites (5.1%)**: 7 peaks

![Differential Composition Charts](09_plots/differential_binding_composition.png)
*Figure 7: Distribution dynamics of high-confidence gained and lost peak selections.*

### 3. Transcriptional Driver Discovery
Motif enrichment calculations reveal distinct regulatory systems steering the cellular response:
- **Genomic Background**: High target alignment to **TEAD3** ($\log_2\text{enr} = 1.64$) and standard AP-1 network complexes (**JUN**, **FOS**, **BATF**).
- **ATAC-seq Background**: Chromatin accessibility-normalized signals isolate structural shifts dominated by **TEAD** family members (TEAD1/2/3/4) and structural anchors (**CTCF**/**CTCFL**).

![Motif Volcano Plot](09_plots/motif_enrichment_volcano.png)
*Figure 5: Volcano dot plot isolating highly significant master regulatory motifs (TEAD3) driving treatment progression.*

### 4. Genomic Relocalization Architecture
Functional peak distribution tracking captures a widespread architectural shift following treatment:
- **Overall Profile**: Dominated by intronic elements (58.8%) and distal intergenic spaces (18.4%).
- **Treatment Impact**: Gained binding centers majorly on intronic regions (61.2%), pointing to distal enhancer activation. Lost binding concentrates inside classic promoters (57.1%), reflecting a selective down-regulation of basal transcription programs.

![Genomic Features Bar Chart](09_plots/genomic_distribution_bar.png)
*Figure 6: Genomic target breakdown comparing localized patterns of gained vs. lost binding domains.*

### 5. Locus Verification Highlights
- **Top Gained Domain (*PDLIM5*)**: Positioned at `chr4:95,393,213–95,393,490`. Displays complete treatment-induced activation, shifting from baseline state (Mean Untreated: 1.0) to strong target capture (Mean Treated: 17.0; $\text{Fold Change} = 9$).
- **Top Lost Domain (*SWI5*/*TRUB2*)**: Positioned at `chr9:131,057,762–131,058,039`. Reflects a dynamic loss of accessibility ($\log_2\text{FC} = -1.58$), closing down an active, baseline euchromatin state upon treatment.

![PDLIM5 Gain Site Verification](09_plots/pdlim5_gain_site.png)
*Figure 10: High-resolution tracking validating structural signal activation across the PDLIM5 locus.*

---

## 🛠️ Computational Requirements
The pipeline was verified on Windows using **R version 4.5.2**.

### Main Bioconductor Infrastructure
```R
library(QuasR)          # Structural Alignment & Framework Management
library(Rsamtools)      # BAM/SAM Quality Management
library(MACSr)          # Interface for peak calling via MACS3
library(monaLisa)       # Binned Motif Enrichment Analytics
library(ChIPseeker)     # Genomic Feature Annotations
```

---

## 📥 Local Workspace Deployment Instructions
To run this pipeline without hitting size limit errors, ensure your project workspace is organized using relative structural references:

1. Clone this project repository into your local terminal directory:
   ```bash
   git clone https://github.com
   ```
2. Open the accompanying `.Rproj` file using RStudio.
3. Manually download your local raw BAM elements into your workspace root folder. The underlying `.gitignore` rules will safeguard your workflow, meaning large data footprints are kept strictly local while script versioning remains lightweight on GitHub.
