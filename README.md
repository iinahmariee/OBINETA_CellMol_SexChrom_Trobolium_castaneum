# Cell and Molecular Biology
## Sex-Chromosome Dosage in a Selected Animal Nuclear Genome Using NCBI, Galaxy, and R

## Objectives

Select one unique animal species; find an NCBI chromosome-level nuclear genome and whole-genome DNA reads; map reads in Galaxy; calculate sex-chromosome-to-autosome depth ratios in R; compare the inference with sample metadata; and document all work on GitHub.

## Part A. Choose and register your animal

| Field | Value |
|---|---|
| Full name | Obiñeta, Inah Marie |
| Section | B |
| Animal species | *Tribolium castaneum* (red flour beetle) |
| Assembly accession | GCF_000002335.3 (Tcas5.) |
| Assembly BioSample | SAMN30135407 |
| Sex system (XY/ZW) | XY |
| WGS run accession | SRR23651710 |
| WGS BioSample | SAMN33549348 |
| Same donor? | No |

## Part B. Retrieve data and map reads in Galaxy

<img width="1918" height="850" alt="01_Galaxy_Fasta SRR_Upload" src="https://github.com/user-attachments/assets/4c088c4f-b6a5-4d04-abe9-7f4e91126013" />

**Figure 1.** Galaxy history showing the uploaded Tcas5.2 reference FASTA and the queued SRA download (fasterq-dump) of the paired-end WGS run SRR23651710.

<img width="1918" height="852" alt="02_Galaxy_Completed_FastQC_forward reverse" src="https://github.com/user-attachments/assets/d1c59cef-6696-4933-8480-b8592855a83c" />

**Figure 2.** Galaxy history showing completed FastQC jobs (Webpage and RawData outputs) for the forward and reverse reads of SRR23651710.

<img width="1912" height="857" alt="03_FastQC_Adapter_Content_Forward" src="https://github.com/user-attachments/assets/7ed6726b-97f4-43ae-9e51-294194d5f738" />

**Figure 3.** FastQC adapter content for the forward reads of SRR23651710 (FastQC on Tribolium castaneum: forward). The adapter content check passed. No Illumina Universal, Illumina Small RNA, Nextera, poly-A or poly-G sequences were detected at any position (about 0%), so no adapter trimming was needed.

<img width="1915" height="862" alt="04_FastQC_Adapter_Content_reverse" src="https://github.com/user-attachments/assets/6af60553-3c02-41d4-8762-5ba795e250f6" />

**Figure 4.** FastQC adapter content for the reverse reads of SRR23651710 (FastQC on Tribolium castaneum: reverse). The adapter content check passed. No adapter or poly-A/poly-G contamination was detected along the reads (about 0% at every position), so no adapter trimming was needed.

<img width="1918" height="853" alt="05_BWA-MEM_BAM_Output" src="https://github.com/user-attachments/assets/ddd40f6f-e238-4054-8576-77be4cc9b845" />

**Figure 5.** BAM output (dataset 21, 598.7 MB) from the Map with BWA-MEM tool (version 0.7.19+galaxy1), aligning SRR23651710 to the icTriCast1.1 reference from the history, sorted by chromosomal coordinates.

<img width="1918" height="852" alt="06_Samtools_Coverage_Output" src="https://github.com/user-attachments/assets/4ffcbdc6-4bc2-415e-a0e5-682d26c2207c" />

**Figure 6.** Samtools coverage table (tabular output, MAPQ ≥ 20, BASEQ ≥ 20) showing mean depth per chromosome and scaffold for SRR23651710 mapped to icTriCast1.1.

## Part C. Label chromosomes for the R analysis

| Accession Number | Chromosome Group |
|---|---|
| NC_004353.4 | Autosomes (A) |
| NT_033777.3 | Autosomes (A) |
| NT_033778.4 | Autosomes (A) |
| NT_033779.5 | Autosomes (A) |
| NT_037436.4 | Autosomes (A) |
| NC_004354.4 | X Chromosome |
| NC_024512.1 | Y Chromosome |

The five autosomal accessions — NC_004353.4, NT_033777.3, NT_033778.4, NT_033779.5, and NT_037436.4 — are grouped together as autosomes because they serve collectively as the reference set for calculating relative sequencing depth and form the consistent autosomal baseline used for comparison against the sex chromosomes. The X chromosome, with accession NC_004354.4, is designated separately so that its coverage can be measured relative to the autosomal average to determine the X:A ratio and evaluate X-chromosome representation in the sample. Similarly, the Y chromosome, identified by accession NC_024512.1, is listed independently to allow calculation of the Y:A ratio and to assess how completely and accurately the Y-chromosome sequence is represented in the sequencing data.

## Part D. Run the R analysis

<img width="1918" height="1018" alt="07_Folder_Galaxy_Output" src="https://github.com/user-attachments/assets/c6ac9de8-8e34-4be7-897a-b36b9a217bcb" />

**Figure 7.** Project folder for the R analysis in Part D. The Tribolium castaneum folder holds the Galaxy output chromosome_coverage.tsv, the chromosome labels chromosome_groups.csv, the script sex_depth_analysis.R, and the results folder the script created.

### R SCRIPT: sex_depth_analysis.R

#Sex-chromosome dosage from whole-genome DNA read depth
#Input: Samtools coverage TSV + student chromosome labels
cov <- read.delim("chromosome_coverage.tsv", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
names(cov)[1] <- "rname"  # Accepts the '#rname' header
if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")

key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)
if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")
key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))
if (anyDuplicated(key$rname)) stop("Repeated rname in key file")
if (!all(key$group %in% c("A", "X", "Y", "Z", "W")))
  stop("Use only A, X, Y, Z or W in the group column")

x <- merge(key, cov[, c("rname", "meandepth")],
           by = "rname", all.x = TRUE, sort = FALSE)
x$meandepth <- suppressWarnings(as.numeric(x$meandepth))
if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")
if ("X" %in% x$group && "Z" %in% x$group)
  stop("Choose either an XY or a ZW species, not both")

if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
  target <- "X"; one_copy <- "XY-like"; two_copies <- "XX-like"
} else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
  target <- "Z"; one_copy <- "ZW-like"; two_copies <- "ZZ-like"
} else stop("Label exactly one X or exactly one Z chromosome")

A <- median(x$meandepth[x$group == "A"])
if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")
x$relative_depth <- x$meandepth / A
ratio <- x$relative_depth[x$group == target]
call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65) one_copy else "Ambiguous"
if (A < 1) call <- "Inconclusive: very low autosomal depth"

#Additional Y or W evidence can support, but not prove, a call
partner <- if (target == "X") "Y" else "W"
partner_ratio <- if (partner %in% x$group)
  x$relative_depth[x$group == partner][1] else NA_real_

#Export results and a simple, reproducible figure
dir.create("results", showWarnings = FALSE)
write.csv(x, "results/chromosome_relative_depth.csv", row.names = FALSE)
png("results/chromosome_depth.png", width = 1000, height = 580)
par(mar = c(8, 4, 3, 1))
barplot(x$relative_depth, names.arg = x$rname, las = 2,
        ylab = "Depth relative to median autosome",
        main = "Autosomal and sex-chromosome depth", cex.names = 0.75)
abline(h = c(0.5, 1.0), lty = 2)
dev.off()

summary <- c(paste("Median autosome depth:", round(A, 3)),
             paste(target, "/ A ratio:", round(ratio, 3)),
             paste(partner, "/ A ratio:", round(partner_ratio, 3)),
             paste("Dosage interpretation:", call),
             "Compare with BioSample sex and discuss limitations.")
writeLines(summary, "results/sex_inference.txt")
cat(paste(summary, collapse = "\n"), "\n")


<img width="1918" height="1020" alt="08_R_Analysis_Output" src="https://github.com/user-attachments/assets/4bc7c835-6c17-43e9-950b-2809d56e6c65" />

**Figure 8.** R analysis of sex-chromosome dosage in Tribolium castaneum (Tcas5.2). 

RStudio run of sex_depth_analysis.R on chromosome_coverage.tsv and chromosome_groups.csv (5 autosomes, the X and one Y scaffold). The median autosome depth was 2.838×. LGX (NC_007416.3) had an X/A ratio of 0.70, so the script’s call was “Ambiguous”, between the XY-like (≤0.65) and XX-like (≥0.80) cutoffs. The Y/A ratio was 0.157, calculated from a single unlocalized Y scaffold (NW_015450210.1, the largest of the 26 LGY scaffolds), since Tcas5.2 has no assembled Y chromosome. This Y value is supporting evidence only.

<img width="1000" height="580" alt="chromosome_depth" src="https://github.com/user-attachments/assets/ece53361-acdb-4484-a7e4-cfdc8da12c8a" />

**Figure 9.** Depth of each chromosome relative to the median autosome.

Bars show mean depth for five autosomes (NC_007418.3, NC_007422.5, NC_007424.3, NC_007420.3, NC_007417.3), the X chromosome LGX (NC_007416.3), and one unlocalized Y scaffold (NW_015450210.1), each divided by the median autosome depth (2.838×). Dashed lines mark 0.5 (one copy) and 1.0 (two copies). The autosomes cluster near 1.0. The X bar sits at about 0.70, between the one-copy and two-copy expectations, so the script’s call is “Ambiguous”. The Y scaffold bar is about 0.16, but it represents only one of 26 small, unlocalized Y fragments, so it is weak supporting evidence only.

## Part E. Interpret your results

<img width="1308" height="900" alt="09_Record_of_WGSsample" src="https://github.com/user-attachments/assets/254806a0-1f2c-40a2-a039-bb71fb1a3369" />

**Figure 10.** NCBI BioSample record for the WGS sample SAMN33549348.

The record describes a Tribolium castaneum (red flour beetle) sample, strain GA2, adult, whole organism (SRA: SRS16920646; BioProject PRJNA606031). The sex attribute is listed as “pooled male and female”. This is the metadata used to compare with the X:A ratio in Part E;

Where the X:A ratio is 0.70, from a median autosome depth of 2.838× and an LGX depth of 1.986×. This falls between the XY-like cutoff (≤0.65) and the XX-like cutoff (≥0.80), so the script's call is "Ambiguous". The WGS BioSample lists the sex as "pooled male and female", so the reads came from a mix of both sexes and the ratio is expected to fall between 0.5 and 1.0. The result is therefore inconclusive for a single-sex call and consistent with a mixed-sex pool. A ratio of 0.70 would fit roughly 40% females, but this is only a rough estimate at such low coverage.

## Questions for the report

**1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?**

I chose the red flour beetle, Tribolium castaneum. It has an XX/XY sex-determination system: females are XX and males are XY, and males carry one X (hemizygous) plus a small, degenerate Y chromosome The Tcas5.2 reference assembly (GCF_000002335.3) supports this. Its assembly report labels one chromosome as the X (LGX, NC_007416.3), lists the Y only as unlocalized scaffolds (LGY), and has nine numbered autosomes (linkage groups 2–10).

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**

The assembly is Tcas5.2 (RefSeq GCF_000002335.3, GenBank GCA_000002335.3), with assembly BioSample SAMN30135407. The WGS run is SRR23651710 (experiment SRX19533913, BioProject PRJNA606031), from BioSample SAMN33549348. The donors are not the same, because the assembly and the reads come from different BioSample records. My result therefore describes the sequenced beetles (a pool of adult GA2-strain beetles), not the beetle used to build the assembly.

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

I used five autosomes: NC_007418.3, NC_007422.5, NC_007424.3, NC_007420.3 and NC_007417.3 (the largest of the linkage-group chromosomes). The X is NC_007416.3 (LGX). As optional supporting evidence, I used NW_015450210.1 as a Y representative, the largest of the 26 unlocalized LGY scaffolds. I verified all of these in the NCBI Tcas5.2 assembly report, which matches each accession to its linkage group (LG2–LG10, LGX, LGY). I excluded the mitochondrial genome and all other unplaced scaffolds. Each name in chromosome_groups.csv matches the rname column of the coverage file exactly.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

The median autosomal depth is 2.838×. The X depth is 1.986×, so the X:A ratio is 0.70. This falls between the XY-like cutoff (0.65 or lower) and the XX-like cutoff (0.80 or higher), so the script’s call is “Ambiguous” and the sex-chromosome complement is inconclusive. The Y scaffold had a Y/A ratio of 0.157, which is only weak supporting evidence, since it comes from one small fragment of an unassembled Y.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**

The BioSample records the sex as “pooled male and female”, so there is no single sex to match. A mix of males (X:A near 0.5) and females (near 1.0) should give a value in between, and 0.70 fits that. It suggests roughly 40% females, but this is only a rough estimate, because coverage is low (about 2.8×) and depth varies between chromosomes. The metadata supports the ambiguous result but does not replace the sequence evidence. Other sources of uncertainty are a repeat-rich X, mapping bias, pseudoautosomal regions, and the fact that the reads come from a different donor than the assembly.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

An assembly is a map of the chromosomes. It shows what chromosomes exist and how long they are, but not how many copies a particular animal carry. That only shows up when reads from the sample are mapped back to the assembly and the X depth is compared with the autosomal depth. Without reads there is nothing to measure. Chromosome lengths or sequences alone can’t give a dosage ratio, and the result also depends on the sample’s own sex and donor, which the assembly doesn’t record.

## Limitations

The main limitation is low coverage. The median autosome depth is only about 2.8×, so the depth of each chromosome is noisy and the X:A ratio of 0.70 is imprecise. The WGS BioSample is also a "pooled male and female" sample, so a single XX or XY ratio was not expected, and the value falls between the XY (about 0.5) and XX (about 1.0) expectations. The reads come from a different BioSample than the reference assembly, so strain or individual differences can affect mapping. Tcas5.2 has no assembled Y chromosome, so the Y is represented only by small unlocalized scaffolds. The Y/A ratio of 0.157 comes from a single scaffold (NW_015450210.1) and is weak supporting evidence only. The X chromosome is repeat-rich, and reads that map poorly or to several places can distort depth, while pseudoautosomal regions shared between the X and Y can also change the ratio. Finally, the X:A ratio sits between the script's cutoffs, which are teaching guidelines and not diagnostic thresholds, so the sex-chromosome complement is inconclusive.

## Methods and resources

- NCBI assembly record: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002335.3/
- NCBI SRA record: https://www.ncbi.nlm.nih.gov/sra/SRX19533913
- NCBI BioSample record: https://www.ncbi.nlm.nih.gov/biosample/SAMN33549348
- NCBI BioProject record: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA606031
- Galaxy Link: https://usegalaxy.org/u/inah.obineta2/h/obineta-sex-dosage-tribolium-castaneum
