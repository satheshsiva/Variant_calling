## Functional Consequence & Population Frequency Annotation
### Annotation with Ensembl VEP
ClinVar only covers variants previously reported in clinical studies. It tells about novel variants, protein-coding impact (missense, frameshift, splice-site), or how common a variant is in healthy populations (gnomAD allele frequency).

For that, we run Ensembl VEP (Variant Effect Predictor). To avoid downloading a 15+ GB offline cache, use VEP's --database flag to query Ensembl's public GRCh37 database over the network (which runs in ~1–2 minutes on Chromosome 20 PASS variants):

1. Install Ensembl VEP via Conda/Mamba (if not already installed):

```Bash
mamba install -c bioconda ensembl-vep
```

2. Run VEP with protein consequence and gnomAD population frequency flags:

```Bash
vep \
  -i pass_variants.vcf.gz \
  -o vep_annotated.vcf.gz \
  --vcf \
  --compress_output bgzip \
  --database \
  --assembly GRCh37 \
  --port 3337 \
  --symbol \
  --canonical \
  --check_existing \
  --sift b \
  --polyphen b \
  --force_overwrite

bcftools index vep_annotated.vcf.gz
```

What those flags add inside the INFO/CSQ tag:
--database: Uses Ensembl variation/consequence databases  
--assembly GRCh37: Specifies reference assembly  
--symbol: Adds gene symbols  
--canonical: Identifies canonical transcripts  
--sift b: Requests SIFT predictions  
--polyphen b: Requests PolyPhen predictions  
--af_gnomade: Requests gnomAD exome frequencies, where supported by the selected VEP data source.  

3. Inspect high-impact or deleterious coding variants using bcftools +split-vep:bcftools has a built-in plugin (+split-vep) that unpacks VEP's CSQ field into readable columns:

```Bash
bcftools +split-vep \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%gnomADe_AF\t%SIFT\t%PolyPhen\n' \
  -i 'IMPACT="HIGH" || IMPACT="MODERATE"' \
  vep_annotated.vcf.gz | head -n 15
```

*Verification*: We should see a clean table listing gene symbols, functional consequences (such as missense_variant or stop_gained), impact severity (MODERATE or HIGH), gnomAD population frequency, and SIFT/PolyPhen scores.

## Annotation with SnpEff
Often, the VEP freezes, meaning VEP stalled right when it started querying variants. The remote MySQL connection to Ensembl's UK server dropped silently (--database mode has no auto-reconnect). Instead of fighting network timeouts with remote VEP, install SnpEff. Its GRCh37 human database is ~120 MB (compared to VEP's 15+ GB cache), runs 100% offline on your laptop, and annotates your entire Chromosome 20 VCF in 10-15 seconds. it is an alternate annotation, it will not completely replace the VEP.

```Bash
# Install SnpEff via conda/mamba
mamba install -c bioconda snpeff -y

# Download the GRCh37.75 database (~120 MB) and annotate all of chr20 in ~15 seconds
snpEff -Xmx4g -v GRCh37.75 pass_variants.vcf.gz | bgzip > snpeff_annotated.vcf.gz
bcftools index snpeff_annotated.vcf.gz
```

*Inspect the functional consequences* (ANN field):
SnpEff writes gene symbols, transcript IDs, protein changes, and impact severity (HIGH, MODERATE, LOW, MODIFIER) into the INFO/ANN tag. Extract HIGH and MODERATE impact variants with bcftools query:

```Bash
bcftools query \
  -i 'INFO/ANN ~ "HIGH" || INFO/ANN ~ "MODERATE"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
  snpeff_annotated.vcf.gz | \
  awk -F'\t' '{split($5, a, "|"); print $1"\t"$2"\t"$3"\t"$4"\t"a[4]"\t"a[2]"\t"a[3]"\t"a[11]}' | \
  head -n 6
```
*Verification*: Within ~15 seconds of running snpEff, the prompt returns and the bcftools query command prints a clean table showing chromosome, position, alleles, Gene Symbol (a[4]), Consequence (a[2], e.g., missense_variant), Impact (a[3]), and Protein Change (a[11], e.g., p.His62Tyr).

```
20      76771   C       T       DEFB125 missense_variant        MODERATE        p.His62Tyr
20      126155  GCAAA   G       DEFB126 frameshift_variant      HIGH    p.Gln55fs
20      126310  ACC     A       DEFB126 frameshift_variant      HIGH    p.Pro106fs
20      168728  T       A       DEFB128 missense_variant        MODERATE        p.Lys27Asn
20      238435  CTGGTCT C       DEFB132 disruptive_inframe_deletion     MODERATE        p.Val7_Leu8del
20      278687  CCGG    C       ZCCHC3  disruptive_inframe_deletion     MODERATE        p.Ala159del
```


## Annotation with gnomAD
GnomAD hosts chromosome-split sites-only VCF files on Google Cloud Storage, we can download just the Chromosome 20 gnomAD sites file and transfer its allele frequency (AF) tag directly into vep_annotated.vcf.gz using bcftools annotate.

1. Download the GRCh37 gnomAD Chromosome 20 Sites VCF & IndexEnsure you pull the GRCh37 (r2.1.1) file so the genomic coordinates and chromosome 20 naming:

```Bash
# Download gnomAD v2.1.1 exomes Chromosome 20 sites VCF (~180 MB) and its tabix index
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.20.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.20.vcf.bgz.tbi
```
(Note: If we also want non-coding/intronic frequencies beyond exomes, gnomAD genomes chr20 is available at the same path replacing exomes with genomes, though the exomes file is much smaller and covers all coding variants.)

2. Annotate vep_annotated.vcf.gz with gnomAD Allele FrequenciesUse bcftools annotate to copy the overall allele frequency (INFO/AF), allele count (INFO/AC), and total allele number (INFO/AN) into your VCF's INFO column:

```Bash
bcftools annotate \
  -a gnomad.exomes.r2.1.1.sites.20.vcf.bgz \
  -c INFO/AF,INFO/AC,INFO/AN \
  vep_annotated.vcf.gz \
  -Oz -o vep_gnomad.vcf.gz

bcftools index vep_gnomad.vcf.gz
# Verification: Check the header of your new file with 
bcftools view -h vep_gnomad.vcf.gz | grep '##INFO=<ID=AF,'
```
```
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">
```

3. Query Combined VEP Consequences + gnomAD FrequenciesNow you can combine bcftools +split-vep with %INFO/AF to view each variant's coding impact alongside its population frequency (and filter for rare variants like INFO/AF < 0.01 or novel variants where INFO/AF="."):

```Bash
bcftools +split-vep \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%IMPACT\t%INFO/AF\t%SIFT\t%PolyPhen\n' \
  -i '(IMPACT="HIGH" || IMPACT="MODERATE") && (INFO/AF<0.01 || INFO/AF=".")' \
  vep_gnomad.vcf.gz | head -n 6
```
*Verification*: We should see a tab-separated list of rare (AF less than 0.01 or . for unseen in gnomAD exomes) HIGH and MODERATE impact variants with their gene symbols and SIFT/PolyPhen predictions.

```
20      76771   C       T       DEFB125,DEFB125 missense_variant,non_coding_transcript_exon_variant     MODERATE,MODIFIER       0.00131262      tolerated(0.44),.       benign(0),.
20      1791760 C       G       RP5-968J1.1     splice_donor_variant&non_coding_transcript_variant      HIGH    .       .       .
20      1895950 C       G       SIRPA,SIRPA,SIRPA       missense_variant,missense_variant,missense_variant      MODERATE,MODERATE,MODERATE      .       deleterious(0.01),deleterious(0.01),deleterious(0.01)   benign(0.437),benign(0.437),benign(0.007)
20      1895951 C       T       SIRPA,SIRPA,SIRPA       missense_variant,missense_variant,missense_variant      MODERATE,MODERATE,MODERATE      .       tolerated(0.39),tolerated(0.39),tolerated(0.41) benign(0.013),benign(0.013),benign(0.306)
20      1895952 T       C       SIRPA,SIRPA,SIRPA       missense_variant,missense_variant,missense_variant      MODERATE,MODERATE,MODERATE      .       tolerated(0.24),tolerated(0.24),tolerated(0.24) benign(0.007),benign(0.007),benign(0.001)
20      2465304 T       A       ZNF343,ZNF343,ZNF343    splice_acceptor_variant,splice_acceptor_variant,splice_acceptor_variant&non_coding_transcript_variant   HIGH,HIGH,HIGH  0.00260073      .,.,.   .,.,.
```

### Merge three annotation into one VCF file
Because vep_gnomad.vcf.gz already contains our Ensembl VEP functional predictions (INFO/CSQ) and gnomAD population frequencies (INFO/AF), we can layer ClinVar's clinical annotations directly onto it with one more bcftools annotate pass. None of the ClinVar tags (CLNSIG, CLNDN, GENEINFO, RS) collide with VEP's CSQ or gnomAD's AF, so all three datasets coexist cleanly in the INFO column.

1. Layer ClinVar onto vep_gnomad.vcf.gz  
Run bcftools annotate using the clinvar.vcf.gz file.

```Bash
bcftools annotate \
  -a clinvar.vcf.gz \
  -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,INFO/RS,ID \
  vep_gnomad.vcf.gz \
  -Oz -o master_annotated.vcf.gz

bcftools index master_annotated.vcf.gz
```
*Verification*: Check that all three annotation sets are registered in the VCF header (CSQ from VEP, AF from gnomAD, and CLNSIG from ClinVar):

```Bash
bcftools view -h master_annotated.vcf.gz | grep -E '##INFO=<ID=(CSQ|AF|CLNSIG),'
```
We should see three matching ##INFO header lines printed.

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|CLIN_SIG|SOMATIC|PHENO">
##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in samples">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Aggregate germline classification for this single variant; multiple values are separated by a vertical bar">
```

2. Query All Three Databases Together in One Table  
Use bcftools +split-vep to unpack VEP's CSQ fields alongside gnomAD's %INFO/AF and ClinVar's %INFO/CLNSIG and %INFO/CLNDN:

```Bash
bcftools +split-vep \
  -f '%CHROM\t%POS\trs%INFO/RS\t%REF\t%ALT\t[%GT]\t%SYMBOL\t%Consequence\t%IMPACT\t%INFO/AF\t%SIFT\t%PolyPhen\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
  -i '(IMPACT="HIGH" || IMPACT="MODERATE") || INFO/CLNSIG!="."' \
  master_annotated.vcf.gz | head -n 6
```
*What this filter does*: (IMPACT="HIGH" || IMPACT="MODERATE") || INFO/CLNSIG!=".": Captures every variant that either alters protein structure (HIGH or MODERATE VEP impact) or has an entry in ClinVar, displaying its sample genotype ([%GT]), gnomAD allele frequency (%INFO/AF), computational pathogenicity (%SIFT, %PolyPhen), and clinical classification (%INFO/CLNSIG) in a single row.

*Verification*: We should see a 14-column tab-separated table combining dbSNP rsID, genotype, gene symbol, VEP consequence, gnomAD frequency, SIFT/PolyPhen scores, and ClinVar clinical significance in each line.

```
20      76771   rs.     C       T       0/1     DEFB125,DEFB125 missense_variant,non_coding_transcript_exon_variant     MODERATE,MODIFIER       0.00131262      tolerated(0.44),.       benign(0),.     .       .
20      126155  rs11467497      GCAAA   G       0/1     DEFB126,DEFB126 frameshift_variant,intron_variant&non_coding_transcript_variant HIGH,MODIFIER   0.146189        .,.     .,.     Benign  not_specified
20      126310  rs11467417      ACC     A       0/1     DEFB126,DEFB126 frameshift_variant,non_coding_transcript_exon_variant   HIGH,MODIFIER   0.552222        .,.     .,.     Benign  not_specified
20      168728  rs.     T       A       1/1     DEFB128,AL360078.1      missense_variant,downstream_gene_variant        MODERATE,MODIFIER       0.579433        tolerated(1),.  benign(0),.     .       .
20      238435  rs.     CTGGTCT C       0/1     DEFB132 inframe_deletion        MODERATE        0.590484        .       .       .       .
20      278687  rs.     CCGG    C       1/1     ZCCHC3,AL034548.1       inframe_deletion,intron_variant&non_coding_transcript_variant   MODERATE,MODIFIER       0.908915        .,.     .,.     .       .
```

### Export annotation into TSC/CSV
```Bash
# 1. Write the column header row
printf "CHROM\tPOS\tRS_ID\tCLINVAR_ID\tREF\tALT\tQUAL\tGENOTYPE\tALLELIC_DEPTH\tTOTAL_DEPTH\tGENE\tCONSEQUENCE\tIMPACT\tCANONICAL\tGNOMAD_AF\tSIFT\tPOLYPHEN\tCLINVAR_SIG\tCLINVAR_DISEASE\n" > master_variants.tsv

# 2. Extract all annotated variants
bcftools +split-vep \
  -s worst \
  -f '%CHROM\t%POS\trs%INFO/RS\t%ID\t%REF\t%ALT\t%QUAL\t[%GT]\t[%AD]\t[%DP]\t%SYMBOL\t%Consequence\t%IMPACT\t%CANONICAL\t%INFO/AF\t%SIFT\t%PolyPhen\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
  master_annotated.vcf.gz >> master_variants.tsv

# 3. Optional: Convert TSV to CSV (replaces existing commas inside annotations with semicolons first so columns don't shift)
awk -F'\t' '{for(i=1;i<=NF;i++){gsub(/,/, ";", $i)}; OFS=","; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' master_variants.tsv > master_variants.csv
```

*Verification*: Check the row count and preview the first 5 lines of your exported table.

```Bash
wc -l master_variants.tsv
head -n 5 master_variants.tsv | column -t -s $'\t'
```
We should see the header row followed by aligned variant rows.

```
CHROM  POS    RS_ID  CLINVAR_ID  REF  ALT  QUAL     GENOTYPE  ALLELIC_DEPTH  TOTAL_DEPTH  GENE     CONSEQUENCE            IMPACT    CANONICAL  GNOMAD_AF  SIFT  POLYPHEN  CLINVAR_SIG  CLINVAR_DISEASE
20     61098  rs.    .           C    T    78.4117  0/1       24,9           34           .        intergenic_variant     MODIFIER  .          .          .     .         .            .
20     61138  rs.    .           C    CT   63.1094  0/1       23,8           31           .        intergenic_variant     MODIFIER  .          .          .     .         .            .
20     61795  rs.    .           G    T    151.363  0/1       18,14          32           .        intergenic_variant     MODIFIER  .          .          .     .         .            .
20     63244  rs.    .           A    C    222.306  0/1       18,24          42           DEFB125  upstream_gene_variant  MODIFIER  .          .          .     .         .            .
```

### Filter for Rare, Damaging, or Pathogenic Candidate Clinical Variants
In clinical variant interpretation, we filter out common benign polymorphisms (GNOMAD_AF less than 0.01, i.e., under 1% population frequency, or unseen .) and retain variants that are either reported as Pathogenic/Likely_pathogenic/VUS in ClinVar OR predicted to severely disrupt protein function (HIGH impact, or MODERATE missense variants flagged as deleterious by SIFT or probably_damaging by PolyPhen).

```Bash
# 1. Write the header row
printf "CHROM\tPOS\tRS_ID\tCLINVAR_ID\tREF\tALT\tQUAL\tGENOTYPE\tALLELIC_DEPTH\tTOTAL_DEPTH\tGENE\tCONSEQUENCE\tIMPACT\tCANONICAL\tGNOMAD_AF\tSIFT\tPOLYPHEN\tCLINVAR_SIG\tCLINVAR_DISEASE\n" > candidate_clinical_variants.tsv

# 2. Filter for rare + damaging/pathogenic candidates using -s worst
bcftools +split-vep \
  -s worst \
  -f '%CHROM\t%POS\trs%INFO/RS\t%ID\t%REF\t%ALT\t%QUAL\t[%GT]\t[%AD]\t[%DP]\t%SYMBOL\t%Consequence\t%IMPACT\t%CANONICAL\t%INFO/AF\t%SIFT\t%PolyPhen\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
  -i '(INFO/AF<0.01 || INFO/AF=".") && (INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Uncertain_significance" || IMPACT="HIGH" || (IMPACT="MODERATE" && (SIFT~"deleterious" || PolyPhen~"damaging")))' \
  master_annotated.vcf.gz >> candidate_clinical_variants.tsv
```
*Verification*: Check how many candidate variants survived prioritization and inspect them:

```Bash
wc -l master_variants.tsv candidate_clinical_variants.tsv
head -n 5 candidate_clinical_variants.tsv | column -t -s $'\t'
```

Because NA12878 is a healthy benchmark sample, the count will drop from tens of thousands of raw calls down to a small, reviewable shortlist of rare coding or VUS/candidate calls.

```
CHROM  POS       RS_ID  CLINVAR_ID  REF  ALT  QUAL     GENOTYPE  ALLELIC_DEPTH  TOTAL_DEPTH  GENE         CONSEQUENCE                                            IMPACT    CANONICAL  GNOMAD_AF   SIFT               POLYPHEN       CLINVAR_SIG  CLINVAR_DISEASE
20     1791760   rs.    .           C    G    159.624  0/1       6,12           18           RP5-968J1.1  splice_donor_variant&non_coding_transcript_variant     HIGH      YES        .           .                  .              .            .
20     1895950   rs.    .           C    G    225.417  1/1       0,25           25           SIRPA        missense_variant                                       MODERATE  .          .           deleterious(0.01)  benign(0.437)  .            .
20     2465304   rs.    4448480     T    A    218.406  0/1       20,20          40           ZNF343       splice_acceptor_variant                                HIGH      YES        0.00260073  .                  .              .            .
20     11790885  rs.    .           C    CTT  228.318  1/1       2,28           30           AL080274.1   splice_acceptor_variant&non_coding_transcript_variant  HIGH      YES        .           .                  .              .            .
```

__END__
