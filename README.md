# Variant Calling and Annotation Pipeline
An end-to-end bioinformatics workflow for identifying, filtering, annotating, prioritizing, and validating genetic variants using open-source tools.

Dataset: GIAB NA12878 (HG001)
Reference Genome: GRCh37 / hg19
Analysis Region: Chromosome 20
Environment: Linux / WSL2

### Workflow Overview
| Step | Analysis                         | Tools                                |
|------|----------------------------------|--------------------------------------|
| 1    | Variant Calling                  | SAMtools, BCFtools                   |
| 2    | Normalization & Technical QC     | BCFtools                             |
| 3    | Functional & Clinical Annotation | Ensembl VEP, SnpEff, ClinVar, gnomAD |
| 4    | Candidate Prioritization         | BCFtools, SIFT, PolyPhen             |
| 5    | Read-Level Verification          | SAMtools tview, IGV                  |

### Workflow
```
Aligned BAM + Reference FASTA
              |
              ▼
     1. Variant Calling
       (BCFtools)
              |
              ▼
  2. Normalization & QC
       (BCFtools)
              |
              ▼
  3. Functional & Clinical
         Annotation
   (VEP / ClinVar / gnomAD)
              |
              ▼
  4. Candidate Prioritization
   (AF / Impact / ClinVar)
              |
              ▼
   5. Read-Level Verification
       (IGV / tview)
              |
              ▼
       Final Candidates
```

### Key Features
Generate raw variants from aligned sequencing reads.
Normalize indels and filter low-quality variants.
Benchmark calls against GIAB high-confidence regions.
Annotate variants with functional consequences, clinical significance, and population frequencies.
Prioritize rare and potentially damaging variants.
Visually inspect candidate variants using aligned sequencing reads.

### Tools & Resources
SAMtools
BCFtools
Ensembl VEP
SnpEff
ClinVar
gnomAD
IGV
Genome in a Bottle

### Repository Structure
```
Variant_calling/
│
├── README.md
├── LICENSE
│
├── Step_1_Variant_Calling/
├── Step_2_Normalization_QC/
├── Step_3_Functional_Annotation/
├── Step_4_Candidate_Prioritization/
├── Step_5_Read_Verification/
│
└── Results/
    ├── VCF_Samples/
    ├── Statistics/
    ├── Annotation_Tables/
    └── IGV_Screenshots/
```
### Final Outputs
Raw and quality-filtered VCF files.
Variant QC statistics and benchmarking metrics.
Functionally and clinically annotated VCF.
Prioritized candidate variant table (TSV/CSV).
IGV screenshots and read-level verification results.

### Notes
This project is intended for educational and research purposes. Variant prioritization does not establish clinical pathogenicity or replace clinical interpretation.

### License
Distributed under the MIT License. See LICENSE for details.