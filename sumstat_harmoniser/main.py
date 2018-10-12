#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Version 2
#
# Harmonise GWAS summary statistics against a reference VCF

import os
import sys
import gzip
import argparse
from copy import deepcopy
from subprocess import Popen, PIPE
from collections import OrderedDict, Counter
from lib.SumStatRecord import SumStatRecord
from lib.VCFRecord import VCFRecord

def main():
    """ Implements main logic.
    """

    # Get args
    global args
    args = parse_args()

    # Intitate handles and counters
    header_written = False
    strand_counter = Counter()
    code_counter = Counter()
    if args.hm_sumstats:
        out_handle = open_gzip(args.hm_sumstats, "wb")

    # Process each row in summary statistics
    for counter, ss_rec in enumerate(yield_sum_stat_records(args.sumstats,
                                                            args.in_sep)):

        # If set to only process 1 chrom, skip none matching chroms
        if args.only_chrom and not args.only_chrom == ss_rec.chrom:
            continue

        # Validate summary stat record
        ret_code = ss_rec.validate_ssrec()
        if ret_code:
            ss_rec.hm_code = ret_code
            strand_counter['Invalid variant for harmonisation'] += 1

        # # DEBUG print progress
        # if counter % 1000 == 0:
        #     print(counter + 1)

        #
        # Load and filter VCF records ------------------------------------------
        #

        # Skip rows that have code 14 (fail validation)
        if not ss_rec.hm_code:

            # Get VCF reference variants for this record
            vcf_recs = get_vcf_records(
                        args.vcf.replace("#", ss_rec.chrom),
                        ss_rec.chrom,
                        ss_rec.pos)

            # Extract the VCF record that matches the summary stat record
            vcf_rec, ret_code = exract_matching_record_from_vcf_records(
                ss_rec, vcf_recs)

            # Set return code when vcf_rec was not found
            if ret_code:
                ss_rec.hm_code = ret_code

            # If vcf record was found, extract some required values
            if vcf_rec:
                # Get alt allele
                vcf_alt = vcf_rec.alt_als[0]
                # Set variant information from VCF file
                ss_rec.hm_rsid = vcf_rec.id
                ss_rec.hm_chrom = vcf_rec.chrom
                ss_rec.hm_pos = vcf_rec.pos
                ss_rec.hm_other_al = vcf_rec.ref_al
                ss_rec.hm_effect_al = vcf_alt

        else:
            vcf_rec = None

        #
        # Harmonise variants ---------------------------------------------------
        #

        # Skip if harmonisation code exists (no VCF record exists or code 14)
        if ss_rec.hm_code:
            strand_counter['No VCF record found'] += 1

        # Harmonise palindromic alleles
        elif is_palindromic(ss_rec.other_al, ss_rec.effect_al):

            strand_counter['Palindormic variant'] += 1
            if args.hm_sumstats:
                ss_rec = harmonise_palindromic(ss_rec, vcf_rec)

        # Harmonise opposite strand alleles
        elif compatible_alleles_reverse_strand(ss_rec.other_al,
                                               ss_rec.effect_al,
                                               vcf_rec.ref_al,
                                               vcf_alt):

            strand_counter['Reverse strand variant'] += 1
            if args.hm_sumstats:
                ss_rec = harmonise_reverse_strand(ss_rec, vcf_rec)

        # Harmonise same forward alleles
        elif compatible_alleles_forward_strand(ss_rec.other_al,
                                               ss_rec.effect_al,
                                               vcf_rec.ref_al,
                                               vcf_alt):

            strand_counter['Forward strand variant'] += 1
            if args.hm_sumstats:
                ss_rec = harmonise_forward_strand(ss_rec, vcf_rec)

        # Should never reach this 'else' statement
        else:
            sys.exit("Error: Alleles were not palindromic, opposite strand, or "
                     "same strand!")

        # Add harmonisation code to counter
        code_counter[ss_rec.hm_code] += 1

        #
        # Write ssrec to output ------------------------------------------------
        #

        if args.hm_sumstats:

            # Add harmonised other allele, effect allele, eaf, beta, or to output
            out_row = OrderedDict()
            out_row["hm_varid"] = vcf_rec.hgvs()[0] if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_rsid"] = ss_rec.hm_rsid if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_chrom"] = ss_rec.hm_chrom if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_pos"] = ss_rec.hm_pos if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_other_allele"] = ss_rec.hm_other_al.str() if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_effect_allele"] = ss_rec.hm_effect_al.str() if vcf_rec and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_beta"] = ss_rec.beta if ss_rec.beta and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_OR"] = ss_rec.oddsr if ss_rec.oddsr and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_OR_lowerCI"] = ss_rec.oddsr_lower if ss_rec.oddsr_lower and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_OR_upperCI"] = ss_rec.oddsr_upper if ss_rec.oddsr_upper and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_eaf"] = ss_rec.eaf if ss_rec.eaf and ss_rec.is_harmonised else args.na_rep_out
            out_row["hm_code"] = ss_rec.hm_code
            # Add other data from summary stat file
            for key in ss_rec.data:
                value = ss_rec.data[key] if ss_rec.data[key] else args.na_rep_out
                out_row[key] = str(value)

            # Write header
            if not header_written:
                outline = args.out_sep.join([str(x) for x in out_row.keys()]) + "\n"
                out_handle.write(outline.encode("utf-8"))
                header_written = True

            # Write row
            outline = args.out_sep.join([str(x) for x in out_row.values()]) + "\n"
            out_handle.write(outline.encode("utf-8"))

    # Close output handle
    if args.hm_sumstats:
        out_handle.close()

    # Write strand_count stats to file
    if args.strand_counts:
        with open_gzip(args.strand_counts, rw='wb') as out_h:
            for key in sorted(strand_counter.keys()):
                out_h.write('{0}\t{1}\n'.format(key, strand_counter[key]).encode('utf-8'))

    # Write outcome code stats to file
    code_table = {
        1:  'Palindromic; Infer strand; Forward strand; Correct orientation; Already harmonised',
        2:  'Palindromic; Infer strand; Forward strand; Flipped orientation; Requires harmonisation',
        3:  'Palindromic; Infer strand; Reverse strand; Correct orientation; Already harmonised',
        4:  'Palindromic; Infer strand; Reverse strand; Flipped orientation; Requires harmonisation',
        5:  'Palindromic; Assume forward strand; Correct orientation; Already harmonised',
        6:  'Palindromic; Assume forward strand; Flipped orientation; Requires harmonisation',
        7:  'Palindromic; Assume reverse strand; Correct orientation; Already harmonised',
        8:  'Palindromic; Assume reverse strand; Flipped orientation; Requires harmonisation',
        9:  'Palindromic; Drop palindromic; Will not harmonise',
        10: 'Forward strand; Correct orientation; Already harmonised',
        11: 'Forward strand; Flipped orientation; Requires harmonisation',
        12: 'Reverse strand; Correct orientation; Already harmonised',
        13: 'Reverse strand; Flipped orientation; Requires harmonisation',
        14: 'Required fields are not known; Cannot harmonise',
        15: 'No matching variants in reference VCF; Cannot harmonise',
        16: 'Multiple matching variants in reference VCF (ambiguous); Cannot harmonise',
        17: 'Palindromic; Infer strand; EAF or reference VCF AF not known; Cannot harmonise',
        18: 'Palindromic; Infer strand; EAF < --maf_palin_threshold; Will not harmonise' }
    if args.hm_statfile:
        with open_gzip(args.hm_statfile, 'wb') as out_h:
            out_h.write('hm_code\tcount\tdescription\n'.encode('utf-8'))
            for key in sorted(code_counter.keys()):
                out_h.write('{0}\t{1}\t{2}\n'.format(key,
                                                     code_counter[key],
                                                     code_table[key]).encode('utf-8') )

    print("Done!")

    return 0


def parse_args():
    """ Parse command line args using argparse.
    """
    parser = argparse.ArgumentParser(description="Summary statistc harmoniser")

    # Input file args
    infile_group = parser.add_argument_group(title='Input files')
    infile_group.add_argument('--sumstats', metavar="<file>",
                        help=('GWAS summary statistics file'), type=str,
                        required=True)
    infile_group.add_argument('--vcf', metavar="<file>",
                        help=('Reference VCF file. Use # as chromosome wildcard.'), type=str, required=True)

    # Output file args
    outfile_group = parser.add_argument_group(title='Output files')
    outfile_group.add_argument('--hm_sumstats', metavar="<file>",
        help=("Harmonised sumstat output file (use 'gz' extension to gzip)"), type=str)
    outfile_group.add_argument('--hm_statfile', metavar="<file>",
        help=("Statistics from harmonisation process output file. Should only be used in conjunction with --hm_sumstats."), type=str)
    outfile_group.add_argument('--strand_counts', metavar="<file>",
        help=("Output file showing number of variants that are forward/reverse/palindromic"), type=str)

    # Harmonisation mode
    mode_group = parser.add_argument_group(title='Harmonisation mode')
    mode_group.add_argument('--palin_mode', metavar="[infer|forward|reverse|drop]",
                        help=('Mode to use for palindromic variants:\n'
                              '(a) infer strand from effect-allele freq, '
                              '(b) assume forward strand, '
                              '(c) assume reverse strand, '
                              '(d) drop palindromic variants'),
                        choices=['infer', 'forward', 'reverse', 'drop'],
                        type=str)

    # Infer strand specific options
    infer_group = parser.add_argument_group(title='Strand inference options',
        description='Options that are specific to strand inference (--palin_mode infer)')
    infer_group.add_argument('--af_vcf_field', metavar="<str>",
                        help=('VCF info field containing alt allele freq (default: AF_NFE)'),
                        type=str, default="AF_NFE")
    infer_group.add_argument('--infer_maf_threshold', metavar="<float>",
                        help=('Max MAF that will be used to infer palindrome strand (default: 0.42)'),
                        type=float, default=0.42)

    # Global column args
    incols_group = parser.add_argument_group(title='Input column names')
    incols_group.add_argument('--chrom_col', metavar="<str>",
                        help=('Chromosome column'), type=str, required=True)
    incols_group.add_argument('--pos_col', metavar="<str>",
                        help=('Position column'), type=str, required=True)
    incols_group.add_argument('--effAl_col', metavar="<str>",
                        help=('Effect allele column'), type=str, required=True)
    incols_group.add_argument('--otherAl_col', metavar="<str>",
                        help=('Other allele column'), type=str, required=True)
    incols_group.add_argument('--beta_col', metavar="<str>",
                        help=('beta column'), type=str)
    incols_group.add_argument('--or_col', metavar="<str>",
                        help=('Odds ratio column'), type=str)
    incols_group.add_argument('--or_col_lower', metavar="<str>",
                        help=('Odds ratio lower CI column'), type=str)
    incols_group.add_argument('--or_col_upper', metavar="<str>",
                        help=('Odds ratio upper CI column'), type=str)
    incols_group.add_argument('--eaf_col', metavar="<str>",
                        help=('Effect allele frequency column'), type=str)
    incols_group.add_argument('--rsid_col', metavar="<str>",
                        help=('rsID column in the summary stat file'), type=str)

    # Global other args
    other_group = parser.add_argument_group(title='Other args')
    other_group.add_argument('--only_chrom', metavar="<str>",
                        help=('Only process this chromosome'), type=str)
    other_group.add_argument('--in_sep', metavar="<str>",
                        help=('Input file column separator [tab|space|comma|other] (default: tab)'),
                        type=str, default='tab')
    other_group.add_argument('--out_sep', metavar="<str>",
                        help=('Output file column separator [tab|space|comma|other] (default: tab)'),
                        type=str, default='tab')
    other_group.add_argument('--na_rep_in', metavar="<str>",
                        help=('How NA  are represented in the input file (default: "")'),
                        type=str, default="")
    other_group.add_argument('--na_rep_out', metavar="<str>",
                        help=('How to represent NA values in output (default: "")'),
                        type=str, default="")
    other_group.add_argument('--chrom_map', metavar="<str>",
                        help=('Map summary stat chromosome names, e.g. `--chrom_map 23=X 24=Y`'),
                        type=str, nargs='+')

    # Parse arguments
    args = parser.parse_args()

    # Convert input/output separators
    args.in_sep = convert_arg_separator(args.in_sep)
    args.out_sep = convert_arg_separator(args.out_sep)

    # Assert that at least one of --hm_sumstats, --strand_counts is selected
    assert any([args.hm_sumstats, args.strand_counts]), \
        "Error: at least 1 of --hm_sumstats, --strand_counts must be selected"

    # Assert that --hm_statfile is only ever used in conjunction with --hm_sumstats
    if args.hm_statfile:
        assert args.hm_sumstats, \
        "Error: --hm_statfile must only be used in conjunction with --hm_sumstats"

    # Assert that mode is selected if doing harmonisation
    if args.hm_sumstats:
        assert args.palin_mode, \
        "Error: '--palin_mode' must be used with '--hm_sumstats'"

    # Assert that inference specific options are supplied
    if args.palin_mode == 'infer':
        assert all([args.af_vcf_field, args.infer_maf_threshold, args.eaf_col]), \
            "Error: '--af_vcf_field', '--infer_maf_threshold' and '--eaf_col' must be used with '--palin_mode infer'"

    # Assert that OR_lower and OR_upper are used both present if any
    if any([args.or_col_lower, args.or_col_upper]):
        assert all([args.or_col_lower, args.or_col_upper]), \
        "Error: '--or_col_lower' and '--or_col_upper' must be used together"

    # Parse chrom_map
    if args.chrom_map:
        try:
            chrom_map_d = dict([pair.split('=') for pair in args.chrom_map])
            args.chrom_map = chrom_map_d
        except ValueError:
            assert False, \
            'Error: --chrom_map must be in the format `--chrom_map 23=X 24=Y`'

    return args

def convert_arg_separator(s):
    ''' Converts [tab|space|comma|other] to a variable
    '''
    if s == 'tab':
        return '\t'
    elif s == 'space':
        return ' '
    elif s == 'comma':
        return ','
    else:
        return s

def exract_matching_record_from_vcf_records(ss_rec, vcf_recs):
    ''' Extracts the vcf record that matches the summary stat record.
    Args:
        ss_rec (SumStatRecord): object containing summary statistic record
        vcf_recs (list of VCFRecords): list containing vcf records
    Returns:
        tuple(
            a single VCFRecord or None,
            output code
            )
    '''

    # Discard if there are no records
    if len(vcf_recs) == 0:
        return (None, 15)

    # Remove alt alleles that don't match sum stat alleles
    for i in range(len(vcf_recs)):
        non_matching_alts = find_non_matching_alleles(ss_rec, vcf_recs[i])
        for alt in non_matching_alts:
            vcf_recs[i] = vcf_recs[i].remove_alt_al(alt)

    # Remove records that don't have any matching alt alleles
    vcf_recs = [vcf_rec for vcf_rec in vcf_recs if vcf_rec.n_alts() > 0]

    # If there are multiple matching records, resolve using rsid
    if len(vcf_recs) > 1:
        vcf_recs = [vcf_rec for vcf_rec in vcf_recs if vcf_rec.id == ss_rec.rsid]

    # Discard ss_rec if there are no valid vcf_recs
    if len(vcf_recs) == 0:
        return (None, 15)

    # Discard ss_rec if there are multiple records
    if len(vcf_recs) > 1:
        return (None, 16)

    # Given that there is now 1 record, use that
    vcf_rec = vcf_recs[0]

    # Discard ssrec if there are multiple matching alleles
    if vcf_rec.n_alts() > 1:
        return (None, 16)

    return (vcf_rec, None)

def harmonise_palindromic(ss_rec, vcf_rec):
    ''' Harmonises palindromic variant
    Args:
        ss_rec (SumStatRecord): object containing summary statistic record
        vcf_rec (VCFRecord): matching vcf record
    Returns:
        harmonised ss_rec
    '''

    # Mode: Infer strand mode
    if args.palin_mode == 'infer':

        # Extract allele frequency if argument is provided
        if args.af_vcf_field and args.af_vcf_field in vcf_rec.info:
            vcf_alt_af = float(vcf_rec.info[args.af_vcf_field][0])
        else:
            ss_rec.hm_code = 17
            return ss_rec

        # Discard if either MAF is greater than threshold
        if ss_rec.eaf:
            if ( af_to_maf(ss_rec.eaf) > args.infer_maf_threshold or
                af_to_maf(vcf_alt_af) > args.infer_maf_threshold ):
                ss_rec.hm_code = 18
                return ss_rec
        else:
            ss_rec.hm_code = 17
            return ss_rec

        # If EAF and alt AF are concordant, then alleles are on forward strand
        if afs_concordant(ss_rec.eaf, vcf_alt_af):

            # If alleles flipped orientation
            if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
                ss_rec.flip_beta()
                ss_rec.is_harmonised = True
                ss_rec.hm_code = 2
                return ss_rec
            # If alleles in correct orientation
            else:
                ss_rec.is_harmonised = True
                ss_rec.hm_code = 1
                return ss_rec

        # Else alleles are on the reverse strand
        else:

            # Take reverse complement of ssrec alleles
            ss_rec.revcomp_alleles()
            # If alleles flipped orientation
            if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
                ss_rec.flip_beta()
                ss_rec.is_harmonised = True
                ss_rec.hm_code = 4
                return ss_rec
            # If alleles in correct orientation
            else:
                ss_rec.is_harmonised = True
                ss_rec.hm_code = 3
                return ss_rec

    # Mode: Assume palindromic variants are on the forward strand
    elif args.palin_mode == 'forward':

        # If alleles flipped orientation
        if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
            ss_rec.flip_beta()
            ss_rec.is_harmonised = True
            ss_rec.hm_code = 6
            return ss_rec
        # If alleles in correct orientation
        else:
            ss_rec.is_harmonised = True
            ss_rec.hm_code = 5
            return ss_rec

    # Mode: Assume palindromic variants are on the reverse strand
    elif args.palin_mode == 'reverse':

        # Take reverse complement of ssrec alleles
        ss_rec.revcomp_alleles()
        # If alleles flipped orientation
        if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
            ss_rec.flip_beta()
            ss_rec.is_harmonised = True
            ss_rec.hm_code = 8
            return ss_rec
        # If alleles in correct orientation
        else:
            ss_rec.is_harmonised = True
            ss_rec.hm_code = 7
            return ss_rec

    # Mode: Drop palindromic variants
    elif args.palin_mode == 'drop':
        ss_rec.hm_code = 9
        return ss_rec

def harmonise_reverse_strand(ss_rec, vcf_rec):
    ''' Harmonises reverse strand variant
    Args:
        ss_rec (SumStatRecord): object containing summary statistic record
        vcf_rec (VCFRecord): matching vcf record
    Returns:
        harmonised ss_rec
    '''
    # Take reverse complement of ssrec alleles
    ss_rec.revcomp_alleles()
    # If alleles flipped orientation
    if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
        ss_rec.flip_beta()
        ss_rec.is_harmonised = True
        ss_rec.hm_code = 13
        return ss_rec
    # If alleles in correct orientation
    else:
        ss_rec.is_harmonised = True
        ss_rec.hm_code = 12
        return ss_rec

def harmonise_forward_strand(ss_rec, vcf_rec):
    ''' Harmonises forward strand variant
    Args:
        ss_rec (SumStatRecord): object containing summary statistic record
        vcf_rec (VCFRecord): matching vcf record
    Returns:
        harmonised ss_rec
    '''
    # If alleles flipped orientation
    if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
        ss_rec.flip_beta()
        ss_rec.is_harmonised = True
        ss_rec.hm_code = 11
        return ss_rec
    # If alleles in correct orientation
    else:
        ss_rec.is_harmonised = True
        ss_rec.hm_code = 10
        return ss_rec

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

def af_to_maf(af):
    """ Converts an allele frequency to a minor allele frequency
    Args:
        af (float or str)
    Returns:
        float
    """
    # Sometimes AF == ".", in these cases, set to 0
    try:
        af = float(af)
    except ValueError:
        af = 0.0

    if af <= 0.5:
        return af
    else:
        return 1 - af

def get_vcf_records(in_vcf, chrom, pos):
    """ Uses tabix to query VCF file. Parses info from record.
    Args:
        in_vcf (str): vcf file
        chrom (str): chromosome
        pos (int): base pair position
    Returns:
        list of VCFRecords
    """
    response = list(tabix_query(in_vcf, chrom, pos, pos))
    return [VCFRecord(line) for line in response]

def yield_sum_stat_records(inf, sep):
    """ Load lines from summary stat file and convert to SumStatRecord class.
    Args:
        inf (str): input file
        sep (str): column separator

    Returns:
        SumStatRecord
    """
    for row in parse_sum_stats(inf, sep):

        # Replace chrom with --chrom_map value
        chrom = row[args.chrom_col]
        if args.chrom_map:
            chrom = args.chrom_map.get(chrom, chrom)
        # Make sumstat class instance
        ss_record = SumStatRecord(chrom,
                                  row[args.pos_col],
                                  row[args.otherAl_col],
                                  row[args.effAl_col],
                                  row.get(args.beta_col, None),
                                  row.get(args.or_col, None),
                                  row.get(args.or_col_lower, None),
                                  row.get(args.or_col_upper, None),
                                  row.get(args.eaf_col, None),
                                  row.get(args.rsid_col, None),
                                  row)
        yield ss_record

def parse_sum_stats(inf, sep):
    """ Yields a line at a time from the summary statistics file.
    Args:
        inf (str): input file
        sep (str): column separator

    Returns:
        OrderedDict: {column: value}
    """
    with open_gzip(inf, "rb") as in_handle:
        # Get header
        header = in_handle.readline().decode("utf-8").rstrip().split(sep)
        # Assert that all column arguments are contained in header
        for arg, value in args.__dict__.items():
            if '_col' in arg and value:
                assert value in header, \
                'Error: --{0} {1} not found in input header'.format(arg, value)
        # Iterate over lines
        for line in in_handle:
            values = line.decode("utf-8").rstrip().split(sep)
            # Replace any na_rep_in values with None
            values = [value if value != args.na_rep_in else None
                      for value in values]
            # Check we have the correct number of elements
            assert len(values) == len(header), 'Error: column length ({0}) does not match header length ({1})'.format(len(values), len(header))
            yield OrderedDict(zip(header, values))

def open_gzip(inf, rw="rb"):
    """ Returns handle using gzip if gz file extension.
    """
    if inf.split(".")[-1] == "gz":
        return gzip.open(inf, rw)
    else:
        return open(inf, rw)

def tabix_query(filename, chrom, start, end):
    """Call tabix and generate an array of strings for each line it returns.
       Author: https://github.com/slowkow/pytabix
    """
    query = '{}:{}-{}'.format(chrom, start, end)
    process = Popen(['tabix', '-f', filename, query], stdout=PIPE)
    for line in process.stdout:
        yield [s.decode("utf-8") for s in line.strip().split()]

def str2bool(v):
    """ Parses argpare boolean input
    """
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':

    main()
