
## Visual Inspection & Read-Level Verification
Once we have prioritized a shortlist of candidate variants (candidate_clinical_variants.tsv), automated database annotations are never trusted blindly. The final step in variant assessment is visual read alignment verification-confirming that the candidate variants are supported by genuine sequencing reads rather than alignment artifacts, homopolymer slippage, or pseudogene mismapping.

We can verify candidate variants directly in the terminal using samtools tview or visually in IGV (Integrative Genomics Viewer) on Windows.

1. Pick the top candidate variant coordinates from your shortlist   
```Bash
awk -F'\t' 'BEGIN {OFS="\t"} NR>1 {print $1":"$2, $5">"$6, "GT="$8, "AD="$9, "Gene="$11, "Conseq="$12, "ClinVar="$18}' candidate_clinical_variants.tsv | head -n 5
```
*What to look*: Note the CHROM:POS coordinate (for example, 20:2465304) and the expected genotype (0/1 heterozygous = ~50% of reads carry the alternate base; 1/1 homozygous = ~100% of reads carry the alternate base).

```
20:1791760      C>G     GT=0/1  AD=6,12 Gene=RP5-968J1.1        Conseq=splice_donor_variant&non_coding_transcript_variant       ClinVar=.
20:1895950      C>G     GT=1/1  AD=0,25 Gene=SIRPA      Conseq=missense_variant ClinVar=.
20:2465304      T>A     GT=0/1  AD=20,20        Gene=ZNF343     Conseq=splice_acceptor_variant  ClinVar=.
20:11790885     C>CTT   GT=1/1  AD=2,28 Gene=AL080274.1 Conseq=splice_acceptor_variant&non_coding_transcript_variant    ClinVar=.
20:12901275     G>A     GT=0/1  AD=19,19        Gene=RP5-1069C8.2       Conseq=splice_donor_variant&non_coding_transcript_variant       ClinVar=.
```

2. Inspect the raw BAM read pileup at that exact position using samtools tview:  
Replace 20:2465304 below with interesting coordinates.

```Bash
samtools tview -p 20:2465304 na12878_chr20.bam chr20.fa
```
(Controls inside samtools tview: Press . to toggle base dots vs. actual nucleotides, g to jump to a new 20:POS coordinate, and q to quit back to your terminal.)



3. Evaluate the read evidence against 4 visual red flags:  
Allelic Balance (AD): For a true 0/1 heterozygous variant, roughly 35%-65% of stacked reads at that column should show the alternate base (e.g., 15 reference bases, 14 alternate bases). If only 2 out of 30 reads show the variant, it is likely a sequencing error.

Strand Balance: Confirm the alternate base appears on both forward reads (uppercase letters like A, T) and reverse reads (lowercase letters like a, t). If every single alternate base is on only one strand, it is likely a strand-bias artifact.

Read Ends: Check that the variant sits in the middle of reads rather than clustered exclusively at the very first or last 2-3 bases of read ends (where base-calling quality drops).

Homopolymer Runs / Repetitive Sequence: Look at the reference sequence row at the top of tview. Indels inside long stretches of identical bases (e.g., AAAAAAAAA or TTTTTTTTT) are frequently polymerase slippage artifacts.

Optional - Inspect in Desktop IGV (Windows):  
The files we need to visualize are are already accessible on our location. Please open IGV on Windows:

Select genome Human (hg19) (which matches GRCh37 coordinates).

Click File -> Load from File... and select both na12878_chr20.bam and master_annotated.vcf.gz.

Paste 20:POS into the search box to inspect read coverage, base quality shading, and variant tracks visually.

__END__