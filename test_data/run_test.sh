#!/usr/bin/env bash
#

# Create index of the vcf
bgzip -c reference_chr1_vcf.testdata.vcf > reference_chr1_vcf.testdata.vcf.gz
tabix -p vcf reference_chr1_vcf.testdata.vcf.gz

# Run harmoniser on test data
mkdir -p output
python ../sumstat_harmoniser.py --sumstats sum_stats.testdata.tsv \
  --vcf reference_chr#_vcf.testdata.vcf.gz \
  --out output/testdata.output.tsv \
  --log output/testdata.log.tsv.gz \
  --rsid_col rsID \
  --chrom_col chrom \
  --pos_col pos \
  --effAl_col effect_allele \
  --otherAl_col other_allele \
  --eaf_col eaf \
  --beta_col beta \
  --maf_palin_threshold 0.42 \
  --af_vcf_field AF_NFE \
  --af_vcf_min 0.001 \
  --infer_strand True \
  --infer_palin True

# Check that output is as expected
#TODO
