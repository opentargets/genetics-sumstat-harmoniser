#!/usr/bin/env bash
#

# Create index of the vcf
bgzip -c test_data/reference_chr1_vcf.testdata.vcf > test_data/reference_chr1_vcf.testdata.vcf.gz
tabix -p vcf test_data/reference_chr1_vcf.testdata.vcf.gz

# Run preliminary on test data
mkdir -p output
../bin/sumstat_harmoniser \
  --sumstats test_data/sum_stats.testdata.tsv \
  --vcf test_data/reference_chr#_vcf.testdata.vcf.gz \
  --preliminary output/testdata.preliminary.tsv \
  --rsid_col rsID \
  --chrom_col chrom \
  --pos_col pos \
  --effAl_col effect_allele \
  --otherAl_col other_allele \
  --eaf_col eaf \
  --beta_col beta \
  --maf_palin_threshold 0.42 \
  --af_vcf_field AF_NFE \
  --palin_mode infer

# Run harmoniser on test data
mkdir -p output
../bin/sumstat_harmoniser \
  --sumstats test_data/sum_stats.testdata.tsv \
  --vcf test_data/reference_chr#_vcf.testdata.vcf.gz \
  --out output/testdata.output.tsv \
  --stats output/testdata.stats.tsv \
  --rsid_col rsID \
  --chrom_col chrom \
  --pos_col pos \
  --effAl_col effect_allele \
  --otherAl_col other_allele \
  --eaf_col eaf \
  --beta_col beta \
  --maf_palin_threshold 0.42 \
  --af_vcf_field AF_NFE \
  --palin_mode infer

# Check that output is as expected
#TODO
