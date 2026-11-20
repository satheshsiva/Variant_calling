## Variant Calling

Install necessary tools in linux env.  
Codes are optimized for WSL2.

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y wget curl bzip2 ca-certificates

wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh

source ~/.bashrc
conda --version
mamba --version

mamba create -n var python=3.11
mamba activate var

mamba install -c bioconda samtools bcftools
```

This is a sample, we don't need actual BAM files for this run. So, I am downloading chr20 alone.

```bash
# 1. Stream ONLY chromosome 20 directly from the remote GIAB BAM
samtools view -b \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam \
  20 > na12878_chr20.bam

# Remove the cached remote .bai file that samtools downloads into your working directory
rm -f RMNISTHS_30xdownsample.bam.bai

# Index your new local chr20 BAM
samtools index na12878_chr20.bam

# 2. Download matching GRCh37 Chromosome 20 reference FASTA
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
samtools faidx hs37d5.fa.gz 20 > chr20.fa
rm -f hs37d5.fa.gz hs37d5.fa.gz.fai hs37d5.fa.gz.gzi
samtools faidx chr20.fa

# Code to directly download chr20 from web
wget [https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz) && gunzip -c Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz > chr20.fa && samtools faidx chr20.fa
```
Now we have both indexed referece fasta and BAM file for chr20. Generate the VCF from Your Local na12878_chr20.bam

```bash
# Once na12878_chr20.bam and chr20.fa are ready, run bcftools to generate your raw VCF:
bcftools mpileup -Ou -f chr20.fa -a AD,DP,SP na12878_chr20.bam | \
  bcftools call -mv -Oz -f GQ -o raw_variants.vcf.gz

    # explicitly specify the variant calling model 
    bcftools mpileup -Ou -f chr20.fa -a AD,DP,SP na12878_chr20.bam | \
     bcftools call -mv -Oz -f GQ --ploidy 2 -o raw_variants.vcf.gz


bcftools index raw_variants.vcf.gz

# Run this command to confirm your VCF contains variant records:
bcftools view -H raw_variants.vcf.gz | head -n 5
```
```
20      61098   .       C       T       78.4117 .       DP=34;VDB=0.282457;SGB=-0.670168;RPBZ=0.813245;MQBZ=0;MQSBZ=0;BQBZ=-0.154626;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=5,19,5,5;MQ=60 GT:PL:DP:SP:AD:GQ       0/1:112,0,208:34:9:24,9:110
20      61138   .       CTTTTTTTTTTTT   CTTTTTTTTTTTTT  63.1094 .       INDEL;IDV=8;IMF=0.258065;DP=31;VDB=0.892692;SGB=-0.651104;RPBZ=-0.903258;MQBZ=0;MQSBZ=0;BQBZ=-2.34848;SCBZ=-0.589768;MQ0F=0;AC=1;AN=2;DP4=3,20,3,5;MQ=60 GT:PL:DP:SP:AD:GQ       0/1:96,0,76:31:8:23,8:81
20      61795   .       G       T       151.363 .       DP=32;VDB=0.917906;SGB=-0.686358;RPBZ=2.01478;MQBZ=0;MQSBZ=0;BQBZ=1.06221;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=11,7,9,5;MQ=60    GT:PL:DP:SP:AD:GQ       0/1:184,0,210:32:0:18,14:127
20      63244   .       A       C       222.306 .       DP=43;VDB=0.882013;SGB=-0.692831;RPBZ=1.28431;MQBZ=0;MQSBZ=0;BQBZ=0.0129639;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=9,9,12,12;MQ=60 GT:PL:DP:SP:AD:GQ       0/1:255,0,209:42:0:18,24:127
20      63799   .       C       T       113.845 .       DP=37;VDB=0.190927;SGB=-0.680642;RPBZ=1.6875;MQBZ=0;MQSBZ=0;BQBZ=-0.431394;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=12,13,7,5;MQ=60  GT:PL:DP:SP:AD:GQ       0/1:147,0,243:37:1:25,12:127
```

__END__