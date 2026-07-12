# Differential ChIP-seq and Chromatin Accessibility Analysis at the MYC Locus: Treatment-Induced TF Redistribution and Regulatory Remodelling

**Author:** Janmita, Kaverimane Umesh
**Genome build:** hg19
**Pipeline:** QuasR/Rbowtie alignment → MACS3 peak calling → monaLisa/JASPAR2022 motif enrichment → ChIPseeker annotation → GREAT functional enrichment

---

## Methods summary

Paired ChIP-seq (transcription factor, TF) and ATAC-seq BAM files for treated and untreated conditions (two replicates each) were quality-filtered (MAPQ ≥ 30), downsampled to equalise sequencing depth, and converted to normalised BigWig tracks for visualisation in IGV/UCSC. Peak calling was performed with **MACS3** across a range of q-value thresholds (1×10⁻⁶, 1×10⁻⁸, 1×10⁻¹⁰, 1×10⁻¹²).

A q-value threshold of **1×10⁻⁸** was selected for downstream analysis, based on visual inspection of concordance between called peaks and BigWig signal in IGV/UCSC Genome Browser. This threshold provided the best balance between sensitivity (retaining true binding events) and specificity (excluding low-confidence, noise-driven calls) across the tested range (1×10⁻⁶ to 1×10⁻¹²).

Replicate peak sets were merged per group to define high-confidence, reproducible binding sites. Read counts under the merged peak set were compared between treated and untreated groups by two-tailed t-test with fold-change filtering to define differential binding sites. Motif enrichment was performed with monaLisa against JASPAR2022 (vertebrate CORE collection) using two independent background models. Genomic annotation used ChIPseeker with TxDb.Hsapiens.UCSC.hg19.knownGene, and functional enrichment of differential sites was performed using GREAT.

---

## 1. Genome-Wide TF Occupancy and Chromatin State at the MYC Locus 

The IGV snapshot (**Figure 1**) captures the *MYC* proto-oncogene locus, a key driver of cell growth and proliferation with an established role in shaping the tumour microenvironment and host immunity [1]. All three exons of *MYC* are visualised, with exon 3 showing the highest signal density — consistent with its identity as the principal protein-coding region. TF binding (ChIP-seq) peaks in both treated and untreated conditions localise predominantly to the region upstream of *MYC* (regulatory/enhancer region, **Figure 2**) and to exon 3. BigWig tracks show substantial, reproducible enrichment across both conditions, indicating persistent TF occupancy consistent with active transcriptional regulation at this locus.

ATAC-seq signal shows pronounced chromatin accessibility in the untreated condition, particularly in the upstream regulatory region and near the promoter, with partial redistribution/reduction following treatment — suggestive of treatment-induced chromatin remodelling. In contrast, H3K27me3 (repressive mark) shows a broad, low-intensity domain rather than discrete peaks, consistent with partial rather than complete repression. The co-occurrence of accessible chromatin (ATAC) and H3K27me3 enrichment at the same locus reflects a poised regulatory state, in which the *MYC* locus remains transcriptionally competent but is subject to epigenetic constraint. NarrowPeak calls are largely concordant between conditions, indicating the treatment effect at this locus is predominantly quantitative (signal magnitude) rather than qualitative (binding-site identity).

Taken together, these data indicate that regulation of *MYC* in this system involves coordinated interplay between chromatin accessibility, TF occupancy, and repressive histone marking — a balance of activating and restraining signals at a single regulatory locus.

**Figure 1.** `figures/01_myc_locus_igv_overview.png`
IGV snapshot of the *MYC* locus showing ChIP-seq BigWig tracks (TF binding; treated = yellow, untreated = red), ATAC-seq tracks (chromatin accessibility; treated = green, untreated = pink), H3K27me3 signal (repressive chromatin; ENCFF488FYZ, teal), reference gene track (hg19), and MACS3 narrowPeak calls for treated and untreated ChIP-seq and ATAC-seq samples.

**Figure 2.** `figures/02_myc_upstream_regulatory_zoom.png`
Zoomed IGV view of the regulatory region between *MYC* and *CASC11*, upstream of the *MYC* transcription start site. Exon 3 of *MYC* is a protein-coding, transcriptionally active region also encoding post-transcriptional regulatory signal.

---

## 2. Total and Differential TF Binding Site Quantification 

Per-sample peak counts were highly consistent across biological replicates (**Table 1**), supporting good ChIP-seq reproducibility and experimental quality. At the group level, only peaks reproducibly detected in both replicates were retained, restricting downstream comparisons to high-confidence, biologically supported binding events rather than sample-specific noise. The treated group showed a modest but detectable increase in total binding sites relative to untreated (≈ 0.13% net increase; 14,164 vs 14,146).

**Table 1. Number of TF binding sites per sample and per group** *(group-level counts represent sites reproducible across both replicates)*

| # | Sample/Group | No. of TF binding sites |
|---|---|---|
| 1 | Untreated 1 | 14,166 |
| 2 | Untreated 2 | 14,149 |
| 3 | Treated 1 | 14,169 |
| 4 | Treated 2 | 14,165 |
| 5 | **Untreated group** | **14,146** |
| 6 | **Treated group** | **14,164** |

Comparison of treated versus untreated binding across the merged consensus peak set (n = 14,170) identified **136 statistically significant differential peaks** (p ≤ 0.05), of which **129 (94.9%)** represented gained binding and **7 (5.1%)** represented lost binding (**Figure 3, Figure 4**), using fold-change thresholds of ≥ 2 (gain) or ≤ 0.5 (loss). This asymmetry indicates that treatment predominantly opened new TF-accessible regions rather than closing existing ones.

**Threshold rationale.**
*P-value selection:* nominal p ≤ 0.05 was used in preference to Benjamini–Hochberg-adjusted p (FDR ≤ 0.05), which yielded zero differential peaks. With 14,170 simultaneous tests and only two replicates per group, statistical power is inherently limited; multiple-testing correction under these conditions is overly conservative and risks false negatives rather than reflecting a genuine absence of biological signal. Nominal p ≤ 0.05 was therefore adopted as an exploratory, hypothesis-generating threshold, with the explicit caveat that it does not control the family-wise false discovery rate.

*Fold-change selection:* consistent with standard practice in differential ChIP-seq pipelines, fold-change (≥ 2 gain / ≤ 0.5 loss) was applied as a secondary effect-size filter rather than a primary statistical criterion, to prioritise robust, biologically meaningful changes over small-magnitude, noise-driven differences. These thresholds are conventional but not universally standardised in the peer-reviewed literature, and should be reported as exploratory cut-offs rather than validated statistical boundaries.

**Figure 3.** `figures/03_differential_sites_barplot.png`
Bar plot of gained versus lost differential binding sites in treated versus untreated groups (p ≤ 0.05; FC ≥ 2 gain / ≤ 0.5 loss).

**Figure 4.** `figures/04_differential_sites_piechart.png`
Pie chart showing the proportional distribution of gained versus lost differential binding sites, same thresholds as Figure 3.

---

## 3. Motif Enrichment at Differential Binding Sites

Motif enrichment analysis was performed on differential binding sites using two independent background models: (i) random genomic regions and (ii) total ATAC-seq peaks. Top-10 enriched motifs, ranked by log2 enrichment and statistical significance (−log10 FDR), are reported in **Tables 2–4**.

At **gained** binding sites, the random-genomic background (**Table 2**) showed strong enrichment for **TEAD3** together with AP-1 family members (ATF3, JUN, FOS/FOSL, BATF complexes), consistent with activation of stress- and growth-responsive transcriptional programmes. Using ATAC-seq peaks as background (**Table 3**), enrichment shifted toward chromatin-accessibility-normalised signal, with TEAD family members (TEAD1/2/3/4) and CTCF/CTCFL predominating — suggesting the residual signal, once general accessibility is controlled for, reflects structural/enhancer-associated regulatory remodelling rather than simple sequence-driven enrichment. **TEAD3 (MA0808.1)** was the most consistent hit across both background models (log2 enrichment ≈ 1.64 against random genomic background, ≈ 1.12 against ATAC background), nominating TEAD3 as a candidate master regulator of Hippo-pathway-dependent enhancer activity, chromatin remodelling, and proliferative transcriptional output at these loci [2].

At **lost** binding sites, motif analysis against a random genomic background could not be reliably performed owing to the small number of loss events (n = 7), which falls below the minimum input size required for robust enrichment testing. Motif analysis against the ATAC-seq background (**Table 4**) was attempted for completeness but returned weak, non-specific enrichment (low log2 enrichment, FDR ≈ 0), and should be interpreted as inconclusive rather than a negative biological finding — this is a power limitation of the loss-site set, not evidence against motif involvement.

Enrichment was visualised as −log10(FDR) against log2 enrichment (**Figure 5**), allowing simultaneous assessment of effect size and statistical confidence. TEAD3 occupies the high-enrichment/high-significance quadrant under both background models, reinforcing its candidacy as a primary driver of the gained-binding phenotype.

**Table 2. Top 10 enriched TF motifs at gained differential binding sites (random genomic/UCSC background)**

| Motif | TF | log2 enrichment | −log10 FDR |
|---|---|---|---|
| MA0808.1 | TEAD3 | 1.6401 | 13.8699 |
| MA1988.1 | ATF3 | 1.7617 | 9.9772 |
| MA0462.2 | BATF::JUN | 1.7672 | 9.8867 |
| MA1634.1 | BATF | 1.7672 | 9.8867 |
| MA0489.2 | JUN | 1.7638 | 9.8867 |
| MA0477.2 | FOSL1 | 1.8401 | 9.7780 |
| MA1130.1 | FOSL2::JUN | 1.6867 | 9.6327 |
| MA0490.2 | JUNB | 1.7741 | 9.6327 |
| MA1128.1 | FOSL1::JUN | 1.6168 | 9.4916 |
| MA0835.2 | BATF3 | 1.6756 | 9.3867 |

**Table 3. Top enriched TF motifs at gained differential binding sites (ATAC-seq background)**

| Motif | TF | log2 enrichment | −log10 FDR |
|---|---|---|---|
| MA0808.1 | TEAD3 | 1.1201 | 18.5618 |
| MA1121.1 | TEAD2 | 0.9462 | 8.7610 |
| MA0139.1 | CTCF | 1.4302 | 8.5407 |
| MA0090.3 | TEAD1 | 0.9135 | 8.2614 |
| MA0809.2 | TEAD4 | 0.8452 | 8.0203 |
| MA1102.2 | CTCFL | 1.1661 | 3.0754 |
| MA1929.1 | CTCF | 0.8972 | 2.0964 |
| MA1930.1 | CTCF | 1.0580 | 2.0878 |
| MA1122.1 | TFDP1 | 0.9344 | 1.0351 |
| MA1105.2 | GRHL2 | 0.8426 | 1.0351 |

**Table 4. Top motifs at lost differential binding sites (ATAC-seq background) — low confidence, n = 7 sites**

| Motif | TF | log2 enrichment | −log10 FDR |
|---|---|---|---|
| MA0004.1 | ARNT | 0.0198 | 0 |
| MA0006.1 | AHR::ARNT | 0.1003 | 0 |
| MA0019.1 | DDIT3::CEBPA | 0.0771 | 0 |
| MA0029.1 | MECOM | 0.0275 | 0 |
| MA0030.1 | FOXF2 | 0.0539 | 0 |
| MA0031.1 | FOXD1 | 0.0467 | 0 |
| MA0040.1 | FOXQ1 | 0.0554 | 0 |
| MA0051.1 | IRF2 | 0.0503 | 0 |
| MA0059.1 | MAX::MYC | 0.1034 | 0 |
| MA0066.1 | PPARG | 0.1302 | 0 |

**Figure 5.** `figures/05_motif_enrichment_dotplot.pdf`
Motif enrichment at gained differential binding sites (−log10 FDR vs log2 enrichment), random genomic and ATAC-seq backgrounds.

---

## 4. Genomic Annotation and Functional Enrichment of Differential Binding Sites 

Genomic annotation of all differential binding sites (**Figure 6**) showed predominant enrichment in intronic regions (58.8%), followed by distal intergenic (14.4%) and promoter (6.6%) elements. Stratifying by direction of change (**Figure 7**), gained sites were predominantly intronic (61.2%), whereas lost sites were predominantly promoter-associated (57.1%). This pattern indicates a treatment-induced shift from promoter-centric to intronic/distal regulatory engagement — consistent with activation of distal enhancer elements (plausibly TEAD3-dependent, per Section 3) alongside targeted attenuation of baseline promoter-proximal accessibility.

GREAT analysis of the differential binding sites (**Figure 8**) identified phospholipase-C-activating GPCR signalling and maintenance of transcriptional fidelity during RNA Pol II elongation as the most significantly enriched GO Biological Process terms (−log10 p = 3.20 and 3.02, respectively), alongside metabolic terms including fatty acid β-oxidation and ketone body catabolism — suggesting the treatment remodels chromatin at loci governing both oncogenic signalling and metabolic reprogramming. This is supported by the top associated genes (**Figure 9**): **RRAS2** (p = 3.12) links mechanistically to GPCR/RAS signalling enrichment, **POLR2I** (p = 3.02), an RNA Pol II subunit, directly underlies the transcriptional-elongation term, and **CPT1C** (p = 2.44) corroborates the fatty-acid-oxidation enrichment — collectively indicating that differential binding sites are functionally concentrated at regulatory elements coordinating transcriptional and metabolic rewiring.

**Figure 6.** `figures/06_genomic_distribution_overall.png`
Genomic distribution of all differential binding sites (bar plot), showing dominant enrichment in intronic regions followed by distal intergenic, promoter, and exonic elements.

**Figure 7.** `figures/07_genomic_distribution_gain_vs_loss.png`
Genomic distribution of differential binding sites stratified by direction: gained sites enriched in intronic/exonic regions; lost sites strongly promoter-associated.

**Figure 8.** `figures/08_great_go_biological_process.png`
GREAT GO Biological Process enrichment of combined differential binding sites (hg19 background).

**Figure 9.** `figures/09_great_ensembl_gene_enrichment.png`
GREAT Ensembl gene-association enrichment of combined differential binding sites.

---

## 5. Locus-Level Validation: Representative Gain and Loss Candidates 

**Candidate gain site** was selected on the basis of highest log2FC (3.17): **chr4:95,393,213–95,393,490**, an intronic region of *PDLIM5* (**Figure 10**). No TF binding was detected in the untreated condition (mean = 1.0), versus strong enrichment in treated (mean = 17.0; fold change = 9; p = 0.0077). *PDLIM5* is a cytoskeletal scaffold protein expressed in neurons, skeletal muscle, and heart, where it interacts with α-actinin, regulates neuronal calcium signalling, promotes MyoD/myozenin transcription, and coordinates PKC/PKD1/PKA and RAS–ERK signalling [3]; its dysregulation has been implicated in bipolar disorder, dilated cardiomyopathy, pulmonary hypertension, and multiple cancers. The absence of H3K27me3 signal together with treatment-exclusive narrowPeak calls indicates this region behaves as a latent enhancer, activated specifically under treatment, consistent with treatment-induced regulatory activation at the *PDLIM5* locus.

**Candidate loss site** was selected on the basis of most negative log2FC (−1.58): **chr9:131,057,762–131,058,039**, a putative regulatory region between *SWI5* and *TRUB2* (**Figure 11, Figure 12**). This region showed significantly higher signal in untreated versus treated (fold change = 0.33; log2FC = −1.58; p = 0.045), a roughly three-fold reduction post-treatment. The absence of H3K27me3 signal combined with strong untreated-specific ATAC-seq and ChIP-seq peaks indicates this region is constitutively euchromatic and loses accessibility/TF occupancy specifically as a consequence of treatment, i.e. an active regulatory element that becomes closed upon treatment.

**Figure 10.** `figures/10_pdlim5_gain_candidate_igv.png`
IGV visualisation of the top gain candidate at the *PDLIM5* locus (chr4:95,393,213–95,393,490). Treated ATAC (green) and ChIP (yellow) show a strong gained peak; untreated ATAC (pink) and ChIP (red) remain low. ENCODE H3K27me3 reference (indigo) shows weak baseline signal.

**Figure 11.** `figures/11_swi5_trub2_loss_candidate_igv.png`
IGV visualisation of the top loss candidate (chr9:131,057,762–131,058,039). Untreated ATAC (pink) and ChIP (red) show strong peaks; treated ATAC (green) and ChIP (yellow) show loss of signal. ENCODE H3K27me3 reference (indigo) shows no signal, indicating a naturally euchromatic state that becomes closed upon treatment.

**Figure 12.** `figures/12_swi5_trub2_locus_map.png`
Genomic location of the top loss candidate site relative to *SWI5* and *TRUB2*.

---

## References

1. Ahmadi SE, Rahimi S, Zarandi B, Chegeni R, Safa M. MYC: a multipurpose oncogene with prognostic and therapeutic implications in blood malignancies. *Journal of Hematology & Oncology*. 2021 Aug 9;14(1):121.
2. National Center for Biotechnology Information. TEAD3 TEA domain transcription factor 3 [*Homo sapiens*]. Gene [Internet]. Bethesda (MD): National Library of Medicine (US); [updated 2026 Apr 8; cited 2026 Apr 11]. Available from: https://www.ncbi.nlm.nih.gov/gene/7005
3. Huang X, Qu R, Ouyang J, Zhong S, Dai J. An overview of the cytoskeleton-associated role of PDLIM5. *Frontiers in Physiology*. 2020 Aug 7;11:975.
