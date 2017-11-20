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
  If no alt alleles remaining -> discard (code: ??) and break

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

import os
import sys
import gzip
from subprocess import Popen, PIPE
from collections import OrderedDict

def main():
    """ Implements main logic
    """

    # Get args
    global args
    args = Args_placeholder()

    # Open output handle
    out_hanlde = open_gzip(args.out_harmonised, "w")
    header_written = False

    # Process each row in summary statistics
    for ss_rec in yield_sum_stat_records(args.in_sum_stats, args.in_sep):


        #
        # Load and filter VCF record -------------------------------------------
        #

        # Get VCF reference variant for this record
        vcf_rec = get_vcf_record(
                    args.in_reference_vcf_pattern.replace("#", ss_rec.chrom),
                    ss_rec.chrom,
                    ss_rec.pos)
        # Discard if there are no records
        if not vcf_rec:
            continue # TODO log that no VCF record was found for this variant
        # Remove alt alleles if allele freq < threshold
        vcf_rec = vcf_rec.filter_alts_by_af(args.af_vcf_min,
                                            args.af_vcf_field)
        # Discard ssrec if there are no alts after filtering
        if vcf_rec.n_alts() == 0:
            continue # TODO log that no alts remained after filtering

        #
        # Remove non-matching multialleleic sites ------------------------------
        #

        # Find and remove non matching alts
        non_matching_alts = find_non_matching_alleles(ss_rec, vcf_rec)
        for alt in non_matching_alts:
            vcf_rec = vcf_rec.remove_alt_al(alt)
        # Discard ssrec if there are no matching alleles
        if vcf_rec.n_alts() == 0:
            continue # TODO log that there were no matching alleles
        # Discard ssrec if there are multiple matching alleles
        if vcf_rec.n_alts() > 1:
            continue # TODO log that there were multiple possibles alleles
        # Given that only 1 alt should exist now, extract alt and af
        vcf_alt = vcf_rec.alt_als[0]
        vcf_alt_af = float(vcf_rec.info[args.af_vcf_field][0])

        #
        # Harmonise variants ---------------------------------------------------
        #

        # Harmonise palindromic alleles
        if is_palindromic(ss_rec.other_al, ss_rec.effect_al):
            if args.infer_palin==True and args.infer_strand==True:

                # Discard if either MAF is greater than threshold
                if ss_rec.eaf:
                    if (af_to_maf(ss_rec.eaf) > args.maf_palin_infer_threshold or
                        af_to_maf(vcf_alt_af) > args.maf_palin_infer_threshold):
                        continue # TODO log that MAFs were greater than threshold
                else:
                    continue # TODO log that no eaf in sumstat file
                # Flip if eafs are not concordant
                if not afs_concordant(ss_rec.eaf, vcf_alt_af):
                    ss_rec.flip_beta()
                    ss_rec.revcomp_alleles()
            else:
                continue # TODO log that it was palindromic

        # Harmonise opposite strand alleles
        elif compatible_alleles_reverse_strand(ss_rec.other_al,
                                               ss_rec.effect_al,
                                               vcf_rec.ref_al,
                                               vcf_alt):
            if args.infer_strand:
                # Take reverse complement of ssrec alleles
                ss_rec.revcomp_alleles()
                # Flip if ss effect allele matches vcf ref allele
                if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
                    ss_rec.flip_beta()
            else:
                continue # TODO log that assuming forward strand, therefore alleles ambiguous

        # Harmonise same strand alleles
        elif compatible_alleles_forward_strand(ss_rec.other_al,
                                               ss_rec.effect_al,
                                               vcf_rec.ref_al,
                                               vcf_alt):
            # Flip if the effect allele matches the vcf ref alleles
            if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
                ss_rec.flip_beta()

        else:
            sys.exit("Error: Alleles were not palindromic, opposite strand, or same strand!")

        #
        # Write ssrec to output ------------------------------------------------
        #

        # Add harmonised other allele, effect allele, eaf, beta to output
        out_row = ss_rec.data
        out_row["OtherAl_hm"] = ss_rec.other_al.str()
        out_row["EffectAl_hm"] = ss_rec.effect_al.str()
        if ss_rec.eaf:
            out_row["EAF_hm"] = str(ss_rec.eaf)
        else:
            out_row["EAF_hm"] = "NA"
        out_row["Beta_hm"] = str(ss_rec.beta)

        # Write header
        if not header_written:
            out_hanlde.write(args.out_sep.join(out_row.keys()) + "\n")
            header_written = True

        # Write row
        out_hanlde.write(args.out_sep.join(out_row.values()) + "\n")

        # print ss_rec # Debug

    # Close handle
    out_hanlde.close()

    print "Done!"

    return 0

def afs_concordant(af1, af2):
    """ Checks whether the allele frequencies of two palindromic variants are
        concordant. Concordant if both are either >0.5 or both are <0.5.
    Args:
        af1, af2 (float): Allele frequencies from two datasets
    Returns:
        Bool: True if concordant
    """
    assert isinstance(af1, float) and isinstance(af2, float)
    if (af1 >= 0.5 and af2 >= 0.5) or (af1 < 0.5 and af2 < 0.5):
        return True
    else:
        return False

def is_palindromic(A1, A2):
    """ Checks if two alleles are palindromic.
    Args:
        A1, A2 (Seq): Alleles (i.e. other and effect alleles)
    """
    return A1.str() == A2.revcomp().str()

def find_non_matching_alleles(sumstat_rec, vcf_rec):
    """ For each vcfrec ref-alt pair check whether it matches either the
        forward or reverse complement of the sumstat alleles.
    Args:
        sumstat_rec (SumStatRecord)
        vcf_rec (VCFRecord)
    Returns:
        list of alt alleles to remove
    """
    alts_to_remove = []
    for ref, alt in vcf_rec.yeild_alleles():
        if not compatible_alleles_either_strand(sumstat_rec.other_al,
                                                sumstat_rec.effect_al,
                                                ref,
                                                alt):
            alts_to_remove.append(alt)
    return alts_to_remove

def compatible_alleles_either_strand(A1, A2, B1, B2):
    """ Checks whether alleles are compatible, either on the forward or reverse
        strand
    Args:
        A1, A2 (Seq): Alleles from one source
        B1, B2 (Seq): Alleles from another source
    Returns:
        Boolean
    """
    return (compatible_alleles_forward_strand(A1, A2, B1, B2) or
            compatible_alleles_reverse_strand(A1, A2, B1, B2))

def compatible_alleles_forward_strand(A1, A2, B1, B2):
    """ Checks whether alleles are compatible on the forward strand
    Args:
        A1, A2 (Seq): Alleles from one source
        B1, B2 (Seq): Alleles from another source
    Returns:
        Boolean
    """
    return set([A1.str(), A2.str()]) == set([B1.str(), B2.str()])

def compatible_alleles_reverse_strand(A1, A2, B1, B2):
    """ Checks whether alleles are compatible on the forward strand
    Args:
        A1, A2 (Seq): Alleles from one source
        B1, B2 (Seq): Alleles from another source
    Returns:
        Boolean
    """
    return set([A1.str(), A2.str()]) == set([B1.revcomp().str(), B2.revcomp().str()])

class VCFRecord:
    """ Holds info from a single vcf row
    """
    def __init__(self, row):
        # Parse fields
        self.chrom = str(row[0])
        self.pos = int(row[1])
        self.id = str(row[2])
        self.ref_al = Seq(row[3])
        self.alt_als = [Seq(nucl) for nucl in row[4].split(",")]
        self.qual = row[5]
        self.filter = str(row[6])
        self.info = parse_info_field(row[7])

    def n_alts(self):
        """ Returns the number of alt alleles """
        return len(self.alt_als)

    def remove_alt_al(self, alt_al):
        """ Given an alt allele, will remove from the record. Will also remove
            corresponding info fields.
        """
        n_alts = len(self.alt_als)
        alt_al_index = self.alt_als.index(alt_al)
        # Remove from alt_al list
        del self.alt_als[alt_al_index]
        # Remove from info fields if same length
        for key in self.info:
            if len(self.info[key]) == n_alts:
                del self.info[key][alt_al_index]
        return self

    def filter_alts_by_af(self, maf_threshold, af_field):
        """ Will remove alt alleles if maf < maf_threshold. Also removes
            corresponding info field.
        Args:
            maf_threshold (float): min threshold to be kept
            af_field (str): INFO field containing AF
        Returns:
            VCF_record
        """
        # Find alts to remove
        alts_to_remove = []
        for alt_al, af in zip(self.alt_als, self.info[af_field]):
            maf = af_to_maf(af)
            if maf < maf_threshold:
                alts_to_remove.append(alt_al)
        # Remove alts
        for alt in alts_to_remove:
            self = self.remove_alt_al(alt)
        return self

    def yeild_alleles(self):
        """ Iterates over alt alleles, yielding list of (ref, alt) tuples.
        Returns:
            List of tuples: (ref allele, alt allele)
        """
        for alt in self.alt_als:
            yield (self.ref_al, alt)

def af_to_maf(af):
    """ Converts an allele frequency to a minor allele frequency
    Args:
        af (float or str)
    Returns:
        float
    """
    af = float(af)
    if af <= 0.5:
        return af
    else:
        return 1 - af

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
    assert len(response) <= 1 # The gnomad vcf should have 1 line per position
    if len(response) == 1:
        vcf_rec = VCFRecord(response[0])
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

    def str(self):
        """ Returns nucleotide as a string
        """
        return str(self.seq)

class SumStatRecord:
    """ Class to hold a summary statistic record.
    """
    def __init__(self, rsid, chrom, pos, other_al, effect_al, beta, eaf, data):
        self.rsid = str(rsid)
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.other_al = Seq(other_al)
        self.effect_al = Seq(effect_al)
        self.beta = float(beta)
        self.data = data
        self.flipped = False
        # Effect allele frequency is not required if we assume +ve strand
        if eaf:
            self.eaf = float(eaf)
            assert 0<= self.eaf <= 1
        else:
            self.eaf = None
        # Assert that chromosome is permissible
        permissible = set([str(x) for x in range(1, 23) + ["X", "Y", "MT"]])
        assert set(self.chrom).issubset(permissible)
        # Assert that other and effect alleles are different
        assert self.other_al.str() != self.effect_al.str()

    def revcomp_alleles(self):
        """ Reverse complement both the other and effect alleles.
        """
        self.effect_al = self.effect_al.revcomp()
        self.other_al = self.other_al.revcomp()

    def flip_beta(self):
        """ Flip the beta, alleles and eaf. Set flipped to True.
        Args:
            revcomp (Bool): If true, will take reverse complement in addition
                            to flipping.
        """
        # Flip beta
        self.beta = self.beta * -1
        # Switch alleles
        new_effect = self.other_al
        new_other = self.effect_al
        self.other_al = new_other
        self.effect_al = new_effect
        # Flip eaf
        if self.eaf:
            self.eaf = 1 - self.eaf
        # Set flipped
        self.flipped = True

    def alleles(self):
        """
        Returns:
            Tuple of (other, effect) alleles
        """
        return (self.other_al, self.effect_al)

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
                                  row.get(args.ss_eaf_col, None),
                                  row
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
