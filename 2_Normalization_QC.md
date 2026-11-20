
## Inspect Summary Metrics & Filter Low-Quality Calls
Now, we have raw_variants.vcf.gz, we assess technical call quality before looking at biological meaning. Raw variant callers output many false positives caused by low read coverage, mapping errors, or strand bias.

First, generate baseline summary statistics on our raw calls, then apply standard quality filters to separate high-confidence variants from sequencing noise.

### 1. Generate raw variant statistics:
```Bash
bcftools stats raw_variants.vcf.gz > raw_stats.txt

# View key summary counts (SNPs, indels, multiallelic sites, and Ti/Tv ratio)
grep "^SN\|^TSTV" raw_stats.txt
```
*What to look for*:
SN (Summary Numbers): Total number of records, SNPs, and indels.
TSTV (Transitions / Transversions ratio): For human whole-genome/chromosome calls, a healthy Ts/Tv ratio is ~2.0-2.1 (and ~2.8-3.0 in exomes). If your raw Ts/Tv is significantly lower (e.g., ~1.5), it indicates excess false-positive calls that need filtering.
```bash
SN      0       number of samples:      1
SN      0       number of records:      109246
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 91729
SN      0       number of MNPs: 0
SN      0       number of indels:       17517
SN      0       number of others:       0
SN      0       number of multiallelic sites:   1138
SN      0       number of multiallelic SNP sites:       31
TSTV    0       61707   30053   2.05    61698   30031   2.05
```

### 2. Normalize indels and apply hard quality filters:
Before filtering, left-align and normalize indels (bcftools norm) so complex variants are represented consistently, then filter out low-confidence calls (QUAL < 30, read depth DP < 10, or excessive strand bias SP > 60):

```Bash
# 1. Normalize variants
bcftools norm -f chr20.fa -m -any \
  raw_variants.vcf.gz \
  -Oz -o normalized_variants.vcf.gz

bcftools index normalized_variants.vcf.gz

# 2. Apply quality filters
bcftools filter -s "LowQual" \
  -e 'QUAL<30 || INFO/DP<10 || FMT/GQ<20 || FMT/SP>60' \
  normalized_variants.vcf.gz \
  -Oz -o filtered_variants.vcf.gz

bcftools index filtered_variants.vcf.gz

# 3. Extract PASS variants
bcftools view -f PASS filtered_variants.vcf.gz \
  -Oz -o pass_variants.vcf.gz

bcftools index pass_variants.vcf.gz
```
*How this works*:
bcftools norm -m -any: Splits multiallelic sites into separate biallelic records and left-aligns indels against chr20.fa.
bcftools filter -s "LowQual": Soft-filters variants failing the thresholds by tagging their FILTER column with LowQual while passing good calls as PASS (preserving records so you can inspect what was excluded).

### 3. Extract only PASS variants and check the improved Ts/Tv ratio:
```Bash
bcftools view -f PASS filtered_variants.vcf.gz -Oz -o pass_variants.vcf.gz
bcftools index pass_variants.vcf.gz

bcftools stats pass_variants.vcf.gz | grep "^SN\|^TSTV"
```
*Verification*:
Compare the TSTV line from raw_stats.txt against pass_variants.vcf.gz. You should see the total variant count drop slightly while the ts/tv ratio increases closer to ~2.0–2.1.
```
SN      0       number of samples:      1
SN      0       number of records:      103926
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 87623
SN      0       number of MNPs: 0
SN      0       number of indels:       16303
SN      0       number of others:       0
SN      0       number of multiallelic sites:   0
SN      0       number of multiallelic SNP sites:       0
TSTV    0       59338   28285   2.10    59338   28285   2.10
```

## Benchmarking Accuracy Against the GIAB Truth Set
To measure how accurate the bcftools pipeline is, we compare the pass_variants.vcf.gz against the official NIST Genome in a Bottle (GIAB) v4.2.1 truth set for NA12878 on Chromosome 20, restricted strictly to GIAB's high-confidence regions (areas of the genome where the truth calls are verified and free of repetitive mapping blind spots).

*Why we use it*: Whenever we build or tweak a new variant-calling pipeline (for example, testing bcftools vs. GATK HaplotypeCaller vs. DeepVariant, or changing quality filters), we compare our pass_variants.vcf.gz against the GIAB truth set (using tools like hap.py or bcftools isec) to calculate Precision (how many of your calls are real), Recall / Sensitivity (what percentage of true variants you caught), and F1 Score.

*When we need it*: Only when benchmarking or validating a pipeline's technical accuracy-not when analyzing real patient or research samples where no truth set exists.

```Bash
# Download NA12878 GRCh37 truth VCF (~130 MB) and index
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi

# Download NA12878 GRCh37 high-confidence regions BED file (~13 MB)
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed
```

Slice both our calls and the GIAB truth set to Chromosome 20 high-confidence regions:
To make a fair comparison, both VCFs must be restricted to the exact same Chromosome 20 high-confidence BED intervals and normalized identically:

```Bash
# 1. Extract Chromosome 20 from the GIAB BED file
awk '$1=="20"' \
  HG001_GRCh37_1_22_v4.2.1_benchmark.bed \
  > giab_chr20_highconf.bed

# 2. Extract and normalize GIAB truth variants inside chr20 high-confidence regions
bcftools view -R giab_chr20_highconf.bed HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz | \
  bcftools norm -f chr20.fa -m -any -Oz -o giab_truth_chr20.vcf.gz
bcftools index giab_truth_chr20.vcf.gz

# 3. Restrict your pass_variants.vcf.gz to the exact same high-confidence regions
bcftools view -R giab_chr20_highconf.bed pass_variants.vcf.gz -Oz -o query_chr20_highconf.vcf.gz
bcftools index query_chr20_highconf.vcf.gz
```

3. Intersect our calls with the GIAB truth set using bcftools isec:bcftools isec splits variants into three buckets:

*Tabspace*0000.vcf: False Positives (FP) - called by our pipeline, but absent in GIAB truth.

*Tabspace*0001.vcf: False Negatives (FN) - present in GIAB truth, but missed by our pipeline.

*Tabspace*0002.vcf: True Positives (TP) - present in both our calls and GIAB truth.

```Bash
bcftools isec -p isec_benchmark query_chr20_highconf.vcf.gz giab_truth_chr20.vcf.gz

# Count FP, FN, and TP records and calculate Precision, Recall, and F1 Score
FP=$(bcftools view -H isec_benchmark/0000.vcf | wc -l)
FN=$(bcftools view -H isec_benchmark/0001.vcf | wc -l)
TP=$(bcftools view -H isec_benchmark/0002.vcf | wc -l)

awk -v tp="$TP" -v fp="$FP" -v fn="$FN" 'BEGIN {
  prec = tp / (tp + fp);
  rec  = tp / (tp + fn);
  f1   = 2 * (prec * rec) / (prec + rec);
  printf "True Positives (TP):  %d\nFalse Positives (FP): %d\nFalse Negatives (FN): %d\nPrecision:            %.4f (%.2f%%)\nRecall (Sensitivity): %.4f (%.2f%%)\nF1 Score:             %.4f\n", tp, fp, fn, prec, prec*100, rec, rec*100, f1
}'
```

*Verification*:
We should see a summary printout showing your exact TP, FP, FN counts alongside Precision, Recall, and F1 Score (typically ~96%-99%+ Precision and Recall for SNPs in high-confidence regions on 30x NA12878 data).

```
True Positives (TP):  78239
False Positives (FP): 2683
False Negatives (FN): 1308
Precision:            0.9668 (96.68%)
Recall (Sensitivity): 0.9836 (98.36%)
F1 Score:             0.9751
```

__END__