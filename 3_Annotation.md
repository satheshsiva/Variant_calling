
## Functional & Clinical Annotation with Public Databases

Now, we have high-confidence variants (pass_variants.vcf.gz), we assess their biological impact and population frequency using public databases. Without annotation, a VCF only tells you where a DNA change happened—not which gene it hits, how common it is in healthy humans, or whether it causes disease.

For human variants, the standard public workflow uses Ensembl VEP (Variant Effect Predictor) or SnpEff / SnpSift combined with ClinVar (clinical pathogenicity) and gnomAD (population allele frequencies).

Because downloading full offline VEP/gnomAD caches takes 15+ GB, we can run Ensembl VEP via its public REST API (--database mode) or annotate directly with ClinVar's lightweight Chromosome 20 VCF (~5 MB). Here is how to run ClinVar annotation locally in under 1 minute using bcftools annotate:

1. Download the official GRCh37 ClinVar VCF and index:(Make sure it matches our GRCh37 assembly and chromosome 20)
```Bash
# Download GRCh37 ClinVar VCF (~150 MB) and its tabix index
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
```
2. Annotate our pass_variants.vcf.gz with ClinVar clinical significance & disease names:
Use bcftools annotate to match your variants by exact chromosome, position, reference allele, and alternate allele, copying ClinVar's clinical significance (CLNSIG), disease name (CLNDN), and gene symbol (GENEINFO) into your VCF's INFO field:

```Bash
# check the columns
bcftools view -h clinvar.vcf.gz | grep -E '^##INFO=<ID=(CLNSIG|CLNDN|GENEINFO|RS)'

# Annotate ClinVar
bcftools annotate \
  -a clinvar.vcf.gz \
  -c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,INFO/RS,ID \
  pass_variants.vcf.gz \
  -Oz -o annotated_clinvar.vcf.gz

bcftools index -f annotated_clinvar.vcf.gz
```

*What those flags do*:
-a clinvar.vcf.gz: Specifies the public annotation database VCF.  
-c INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO,ID: Transfers Clinical Significance (e.g., Benign, Pathogenic, Uncertain_significance), Condition/Disease Name, Gene Symbol, and rsID (ID column) wherever an exact variant match occurs.

3. Query and inspect annotated variants in a clean table:
Use bcftools query to extract only variants that matched a ClinVar entry (INFO/CLNSIG!=".") into readable columns:

```Bash
bcftools query \
  -i 'INFO/CLNSIG!="."' \
  -f '%CHROM\t%POS\trs%INFO/RS\t%ID\t%REF\t%ALT\t[%GT]\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
  annotated_clinvar.vcf.gz | head -n 6
```

*Verification*:
We should see a tab-separated output showing rsIDs (e.g., rs6054257), the affected gene symbol, the sample's genotype (0/1 heterozygous or 1/1 homozygous alt), and its ClinVar classification (most variants in NA12878 will be Benign or Likely_benign, as NA12878 is a healthy reference individual).

```
20      126155  rs11467497      402585  GCAAA   G       0/1     DEFB126:81623   Benign  not_specified
20      126310  rs11467417      402586  ACC     A       0/1     DEFB126:81623   Benign  not_specified
20      371972  rs.     4944047 T       C       0/1     TRIB3:57761     Benign  not_specified
20      377226  rs.     4944050 C       T       0/1     TRIB3:57761     Benign  not_specified
20      389456  rs3746792       1185505 C       T       0/1     RBCK1:10616     Benign  not_specified|not_provided|Polyglucosan_body_myopathy_type_1
20      391170  rs7059  2628148 C       T       0/1     RBCK1:10616     Benign/Likely_benign    not_specified
```

__END__