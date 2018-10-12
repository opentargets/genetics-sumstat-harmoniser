#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import subprocess as sp
import hashlib

def main():

    # Args
    in_test_vcf = 'test_data/reference_vcf.testdata.vcf'
    in_test_sumstats = 'test_data/sum_stats.testdata.tsv'
    hm_script = '../bin/sumstat_harmoniser'
    stdout_log = 'output/stdout.txt'
    outdir = 'output'
    check_output = True

    # Make output direcctory
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Open stdout redirect
    stdout_h = open(stdout_log, 'w')

    # Make tabix index of vcf
    prepare_vcf_tabix(in_test_vcf)

    # Test using `--strand_counts` on test data
    out_strand_counts = os.path.join(outdir, 'test.strand_counts.tsv')
    expected_strand_counts = 'expected_output/test.strand_counts.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--strand_counts {out_file} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--otherAl_col other_allele').format(
            script=hm_script,
            in_sumstats=in_test_sumstats,
            in_vcf=in_test_vcf,
            out_file=out_strand_counts)
    # print(cmd)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if md5(out_strand_counts) == md5(expected_strand_counts):
            print('Testing --strand_counts: PASS')
        else:
            print('Testing --strand_counts: FAIL')

    # Test harmonise using `--palin_mode infer` on test data
    out_infer_harm = os.path.join(outdir, 'test.palin_infer.harmonised.tsv')
    out_infer_stats = os.path.join(outdir, 'test.palin_infer.stats.tsv')
    expected_infer_harm = 'expected_output/test.palin_infer.harmonised.tsv'
    expected_infer_stats = 'expected_output/test.palin_infer.stats.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--hm_sumstats {out_harm} '
           '--hm_statfile {out_stats} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--otherAl_col other_allele '
           '--eaf_col eaf '
           '--beta_col beta '
           '--or_col OR '
           '--or_col_lower OR_lowerCI '
           '--or_col_upper OR_upperCI '
           '--infer_maf_threshold 0.42 '
           '--af_vcf_field AF_NFE '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--palin_mode infer').format(
           script=hm_script,
           in_sumstats=in_test_sumstats,
           in_vcf=in_test_vcf,
           out_harm=out_infer_harm,
           out_stats=out_infer_stats)
    # print(cmd)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if ((md5(out_infer_harm) == md5(expected_infer_harm)) and
            (md5(out_infer_stats) == md5(expected_infer_stats))):
            print('Testing --palin_mode infer: PASS')
        else:
            print('Testing --palin_mode infer: FAIL')

    # Test harmonise using `--palin_mode forward` on test data
    out_forward_harm = os.path.join(outdir, 'test.palin_forward.harmonised.tsv')
    out_forward_stats = os.path.join(outdir, 'test.palin_forward.stats.tsv')
    expected_forward_harm = 'expected_output/test.palin_forward.harmonised.tsv'
    expected_forward_stats = 'expected_output/test.palin_forward.stats.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--hm_sumstats {out_harm} '
           '--hm_statfile {out_stats} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--otherAl_col other_allele '
           '--eaf_col eaf '
           '--beta_col beta '
           '--or_col OR '
           '--or_col_lower OR_lowerCI '
           '--or_col_upper OR_upperCI '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--palin_mode forward').format(
           script=hm_script,
           in_sumstats=in_test_sumstats,
           in_vcf=in_test_vcf,
           out_harm=out_forward_harm,
           out_stats=out_forward_stats)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if ((md5(out_forward_harm) == md5(expected_forward_harm)) and
            (md5(out_forward_stats) == md5(expected_forward_stats))):
            print('Testing --palin_mode forward: PASS')
        else:
            print('Testing --palin_mode forward: FAIL')

    # Test harmonise using `--palin_mode reverse` on test data
    out_reverse_harm = os.path.join(outdir, 'test.palin_reverse.harmonised.tsv')
    out_reverse_stats = os.path.join(outdir, 'test.palin_reverse.stats.tsv')
    expected_reverse_harm = 'expected_output/test.palin_reverse.harmonised.tsv'
    expected_reverse_stats = 'expected_output/test.palin_reverse.stats.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--hm_sumstats {out_harm} '
           '--hm_statfile {out_stats} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--otherAl_col other_allele '
           '--eaf_col eaf '
           '--beta_col beta '
           '--or_col OR '
           '--or_col_lower OR_lowerCI '
           '--or_col_upper OR_upperCI '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--palin_mode reverse').format(
           script=hm_script,
           in_sumstats=in_test_sumstats,
           in_vcf=in_test_vcf,
           out_harm=out_reverse_harm,
           out_stats=out_reverse_stats)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if ((md5(out_reverse_harm) == md5(expected_reverse_harm)) and
            (md5(out_reverse_stats) == md5(expected_reverse_stats))):
            print('Testing --palin_mode reverse: PASS')
        else:
            print('Testing --palin_mode reverse: FAIL')

    # Test harmonise using `--palin_mode drop` on test data
    out_drop_harm = os.path.join(outdir, 'test.palin_drop.harmonised.tsv')
    out_drop_stats = os.path.join(outdir, 'test.palin_drop.stats.tsv')
    expected_drop_harm = 'expected_output/test.palin_drop.harmonised.tsv'
    expected_drop_stats = 'expected_output/test.palin_drop.stats.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--hm_sumstats {out_harm} '
           '--hm_statfile {out_stats} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--otherAl_col other_allele '
           '--eaf_col eaf '
           '--beta_col beta '
           '--or_col OR '
           '--or_col_lower OR_lowerCI '
           '--or_col_upper OR_upperCI '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--palin_mode drop').format(
           script=hm_script,
           in_sumstats=in_test_sumstats,
           in_vcf=in_test_vcf,
           out_harm=out_drop_harm,
           out_stats=out_drop_stats)
    # print(cmd)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if ((md5(out_drop_harm) == md5(expected_drop_harm)) and
            (md5(out_drop_stats) == md5(expected_drop_stats))):
            print('Testing --palin_mode drop: PASS')
        else:
            print('Testing --palin_mode drop: FAIL')

    # Test harmonise using `--chrom_map 23=X` with `--palin_mode infer` on test data
    out_map_harm = os.path.join(outdir, 'test.chrom_map.harmonised.tsv')
    out_map_stats = os.path.join(outdir, 'test.chrom_map.stats.tsv')
    expected_map_harm = 'expected_output/test.chrom_map.harmonised.tsv'
    expected_map_stats = 'expected_output/test.chrom_map.stats.tsv'
    cmd = ('{script} '
           ' --sumstats {in_sumstats} '
           '--vcf {in_vcf}.gz '
           '--hm_sumstats {out_harm} '
           '--hm_statfile {out_stats} '
           '--chrom_col chrom '
           '--pos_col pos '
           '--effAl_col effect_allele '
           '--otherAl_col other_allele '
           '--eaf_col eaf '
           '--beta_col beta '
           '--or_col OR '
           '--or_col_lower OR_lowerCI '
           '--or_col_upper OR_upperCI '
           '--infer_maf_threshold 0.42 '
           '--af_vcf_field AF_NFE '
           '--chrom_map 23=X '
           '--na_rep_in NA '
           '--na_rep_out NA '
           '--rsid_col rsID '
           '--palin_mode infer').format(
           script=hm_script,
           in_sumstats=in_test_sumstats,
           in_vcf=in_test_vcf,
           out_harm=out_map_harm,
           out_stats=out_map_stats)
    # print(cmd)
    sp.call(cmd, shell=True, stdout=stdout_h)
    if check_output:
        if ((md5(out_map_harm) == md5(expected_map_harm)) and
            (md5(out_map_stats) == md5(expected_map_stats))):
            print('Testing --chrom_map 23=X: PASS')
        else:
            print('Testing --chrom_map 23=X: FAIL')

    # Close stdout handle
    stdout_h.close()

    return 0

def prepare_vcf_tabix(in_vcf):
    ''' Bgzip and tabix the test vcf file
    '''

    # Bgzip
    vcfgz = '{0}.gz'.format(in_vcf)
    cmd = 'bgzip -c {in_vcf} > {out_vcfgz}'.format(
        in_vcf=in_vcf, out_vcfgz=vcfgz)
    sp.call(cmd, shell=True)

    # Tabix
    cmd = 'tabix -p vcf {in_vcfgz}'.format(in_vcfgz=vcfgz)
    sp.call(cmd, shell=True)

    return 0

def md5(fname):
    ''' Produce a md5 hash of a file. Warning, only works with files approx
        <1Mb. https://stackoverflow.com/a/3431838
    '''
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

if __name__ == '__main__':

    main()
