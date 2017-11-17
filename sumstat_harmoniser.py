#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Version 0.1
#
# Harmonise GWAS summary statistics against a reference VCF

"""
# Pseudo code:

Load summary stat record (ssrec).

Get record from reference VCF (vcfrec) using ssrec chrom:pos and tabix.
  If no record -> discard (code: ??) and break
  Assert that only 1 record was returned (should be the case for gnomad vcfs at least)
  Discard alt alleles if allele freq is below af_vcf_min threshold

As vcfrec could be multialleleic, compare ssrec alleles to vcfrec alleles.
Search for ssrec allele match with each vcfrec ref-alt allele pair and their reverse complements.
Remove alt allele and corresponding INFO entries if not a match.
  If no matches -> discard ssrec as ambiguous (code: ??) and break
  If >1 match -> discard ssrec as ambiguous (code: ??) and break
  If 1 match -> proceed

Check whether ssrec alleles and vcfrec alleles are palindromic
  If palindromic
    If infer_palin==True and infer_strand==True:
      If MAF of either ssrec or vcfrec alleles is > maf_palin_infer_threshold -> discard (error code: ??)
      If ssrec and vcfrec EAFs are not concordant (one is >0.5 and other <0.5)
        Flip beta (code: ??) and break
      Else:
        Don't flip beta (code: ??) and break
    Else: -> discard (code: ??) and break

Check whether on the opposite strand: ssrec alleles == reverse_complement(vcfrec alleles)
  If opposite strand
    If infer_strand==True
      If ssrec effect_al == reverse_complement(vcfrec effect_al):
        Don't flip beta (code: ??) and break
      Else:
        Flip beta (code: ??) and break
    Else:
      Assuming forward strand therefore alleles are ambiguous -> discard (code: ??) and break

Check that alleles are on the same strand (they should be by now).
  If same stand:
    If ssrec effect_al == vcfrec effect_al:
      Don't flip beta (code: ??) and break
    Else:
      Flip beta (code: ??) and break

Else error (shouldn't ever get here):
  Raise error
  
"""

import sys
import os
from subprocess import Popen, PIPE
import gzip
from collections import OrderedDict

def main():

    # Get args
    global args
    args = Args_placeholder()

    # Process each row in summary statistics
    for ss_rec in yield_sum_stat_records(args.in_sum_stats, args.in_sep):
        print ss_rec
        # Get VCF reference variant for this record
        vcf_rec = get_vcf_record(
            args.in_reference_vcf_pattern.replace("#", ss_rec.chrom),
            ss_rec.chrom,
            ss_rec.pos)
        if vcf_rec:
            print vcf_rec.info
    return 0

class VCF_record:
    """ Holds info from a single vcf row
    """
    def __init__(self, row):
        # Parse fields
        self.chrom = str(row[0])
        self.pos = int(row[1])
        self.id = str(row[2])
        self.ref_al = str(row[3])
        self.alt_als = row[4].split(",")
        self.qual = row[5]
        self.filter = str(row[6])
        self.info = parse_info_field(row[7])
        # Derive vars
        self.nalts = len(self.alt_als)

def parse_info_field(field):
    """ Parses a vcf info field in to a dictionary
    Args:
        field (str): info field from vcf
    Returns:
        Dict: {name:[val1, val2]}
    """
    d = {}
    for entry in field.split(";"):
        key, value = entry.split("=")
        d[key] = value.split(",")
    return d

def get_vcf_record(in_vcf, chrom, pos):
    """ Uses tabix to query VCF file. Parses info from record.
    Args:
        in_vcf (str): vcf file
        chrom (str): chromosome
        pos (int): base pair position
    Returns:
        VCF_record
    """
    response = list(tabix_query(in_vcf, chrom, pos, pos))
    # The gnomad vcf should have 1 line per position
    assert len(response) <= 1
    if len(response) == 1:
        vcf_rec = VCF_record(response[0])
        return vcf_rec
    else:
        return None

class Seq:
    """ DNA nucleotide class
    Args:
        seq (str): dna sequence
    """
    def __init__(self, seq):
        self.seq = str(seq).upper()
        assert set(list(self.seq)).issubset(set(["A", "T", "G", "C", "D", "I", "N"]))

    def __repr__(self):
        return str(self.seq)

    def rev(self):
        """ Reverse sequence """
        return Seq(self.seq[::-1])

    def comp(self):
        """ Complementary sequence """
        com_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
        com_seq = "".join([com_dict.get(nucl, "N") for nucl in self.seq])
        return Seq(com_seq)

    def revcomp(self):
        """ Reverse complement sequence """
        return Seq(self.rev().comp())

class SumStatRecord:
    """ Class to hold a summary statistic record.
    """
    def __init__(self, rsid, chrom, pos, other_al, effect_al, beta, eaf):
        self.rsid = str(rsid)
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.other_al = Seq(other_al)
        self.effect_al = Seq(effect_al)
        self.beta = float(beta)
        # Effect allele frequency is not required if we assume +ve strand
        if eaf:
            self.eaf = float(eaf)
            assert 0<= self.eaf <= 1
        else:
            self.eaf = None
        # Assert that chromosome is permissible
        permissible = set([str(x) for x in range(1, 23) + ["X", "Y", "MT"]])
        assert set(self.chrom).issubset(permissible)

    def __repr__(self):
        return "\n".join(["Sum stat record",
                          "  rsID         : " + self.rsid,
                          "  chrom        : " + self.chrom,
                          "  pos          : " + str(self.pos),
                          "  other allele : " + str(self.other_al),
                          "  effect allele: " + str(self.effect_al),
                          "  beta         : " + str(self.beta),
                          "  EAF          : " + str(self.eaf)
                          ])


def yield_sum_stat_records(inf, sep):
    """ Load lines from summary stat file and convert to SumStatRecord class.
    Args:
        inf (str): input file
        sep (str): column separator

    Returns:
        SumStatRecord
    """
    for row in parse_sum_stats(inf, sep):
        ss_record = SumStatRecord(row[args.ss_rsid_col],
                                  row[args.ss_chrom_col],
                                  row[args.ss_pos_col],
                                  row[args.ss_otherAl_col],
                                  row[args.ss_effectAl_col],
                                  row[args.ss_beta_col],
                                  row.get(args.ss_eaf_col, None)
                                  )
        yield ss_record

def parse_sum_stats(inf, sep):
    """ Yields a line at a time from the summary statistics file.
    Args:
        inf (str): input file
        sep (str): column separator

    Returns:
        OrderedDict: {column: value}
    """
    with open_gzip(inf, "r") as in_handle:
        header = in_handle.readline().rstrip().split(sep)
        for line in in_handle:
            values = line.rstrip().split(sep)
            assert(len(values) == len(header))
            yield OrderedDict(zip(header, values))

def open_gzip(inf, rw="r"):
    """ Returns handle using gzip if gz file extension.
    """
    if inf.split(".")[-1] == "gz":
        return gzip.open(inf, rw)
    else:
        return open(inf, rw)

class Args_placeholder:
    """ Placeholder for argparse.
    """
    def __init__(self):
        # File args
        self.name = "testdata"
        self.in_sum_stats = "test_data/sum_stats.testdata.tsv"
        self.in_reference_vcf_pattern = "test_data/reference_chr#_vcf.testdata.vcf.gz"
        self.out_harmonised = "output/output.testdata.tsv"
        self.out_log_prefix = "logs/testdata"
        # Column args
        self.ss_rsid_col = "rsID"
        self.ss_chrom_col = "chrom"
        self.ss_pos_col = "pos"
        self.ss_effectAl_col = "effect_allele"
        self.ss_otherAl_col = "other_allele"
        self.ss_eaf_col = "eaf"
        self.ss_beta_col = "beta"
        # Other args
        self.maf_palin_infer_threshold = 0.42
        self.af_vcf_field = "AF_NFE"
        self.af_vcf_min = 0.001 # 0.1% MAF
        self.in_sep = "\t"
        self.out_sep = "\t"
        self.infer_strand = True
        self.infer_palin = True
        self.infer_multialt = True

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns.
       Author: https://github.com/slowkow/pytabix
    """
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        yield line.strip().split()

if __name__ == '__main__':

    main()
