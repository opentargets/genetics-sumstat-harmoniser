#!/usr/bin/env bash
#

# Create index of the vcf
bgzip -c test_data/reference_chr1_vcf.testdata.vcf > test_data/reference_chr1_vcf.testdata.vcf.gz
tabix -p vcf test_data/reference_chr1_vcf.testdata.vcf.gz

# Test `--preliminary` on test data
mkdir -p output
../bin/sumstat_harmoniser \
  --sumstats test_data/sum_stats.testdata.tsv \
  --vcf test_data/reference_chr#_vcf.testdata.vcf.gz \
  --strand_counts output/testdata.strand_counts.tsv \
  --chrom_col chrom \
  --pos_col pos \
  --effAl_col effect_allele \
  --otherAl_col other_allele

# Run harmoniser on test data
mkdir -p output
../bin/sumstat_harmoniser \
  --sumstats test_data/sum_stats.testdata.tsv \
  --vcf test_data/reference_chr#_vcf.testdata.vcf.gz \
  --hm_sumstats output/testdata.output.tsv \
  --hm_statfile output/testdata.stats.tsv \
  --chrom_col chrom \
  --pos_col pos \
  --effAl_col effect_allele \
  --otherAl_col other_allele \
  --eaf_col eaf \
  --beta_col beta \
  --infer_maf_threshold 0.42 \
  --af_vcf_field AF_NFE \
  --palin_mode infer

# Check that output is as expected
#TODO
