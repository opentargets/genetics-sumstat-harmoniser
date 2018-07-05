#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Version 0.1
#
# Harmonise GWAS summary statistics against a reference VCF

import os
import sys
import gzip
import argparse
from copy import deepcopy
from subprocess import Popen, PIPE
from collections import OrderedDict, Counter

def main():
    """ Implements main logic.
    """

    # Get args
    global args
    args = parse_args()

    # Open output handles
    out_hanlde = open_gzip(args.out, "wb")
    log_hanlde = open_gzip(args.log, "wb")
    header_written = False
    stats = initiate_stats()

    # Process each row in summary statistics
    for counter, ss_rec in enumerate(yield_sum_stat_records(args.sumstats, args.in_sep)):

        # If set to only process 1 chrom, skip none matching chroms
        if args.only_chrom and not args.only_chrom == ss_rec.chrom:
            continue

        # Make a unmodified copy of the summary statistic record
        ss_rec_raw = deepcopy(ss_rec)

        # DEBUG print progress
        if counter % 1000 == 0:
            print(counter + 1)

        #
        # Load and filter VCF records ------------------------------------------
        #

        # Get VCF reference variants for this record
        vcf_recs = get_vcf_records(
                    args.vcf.replace("#", ss_rec.chrom),
                    ss_rec.chrom,
                    ss_rec.pos)

        # Discard if there are no records
        if len(vcf_recs) == 0:
            stats["Pre-filter"]["No record in VCF, discarded"] += 1
            write_to_log(log_hanlde, ss_rec_raw, "None; No record in VCF; Discarded")
            continue

        # Remove alt alleles if allele freq < threshold
        vcf_recs = [vcf_rec.filter_alts_by_af(args.af_vcf_min, args.af_vcf_field)
                    for vcf_rec in vcf_recs]

        # Remove alt alleles that don't match sum stat alleles
        for i in range(len(vcf_recs)):
            non_matching_alts = find_non_matching_alleles(ss_rec, vcf_recs[i])
            for alt in non_matching_alts:
                vcf_recs[i] = vcf_recs[i].remove_alt_al(alt)

        # Remove records that don't have any matching alt alleles
        vcf_recs = [vcf_rec for vcf_rec in vcf_recs if vcf_rec.n_alts() > 0]

        # Discard ss_rec if there are no valid vcf_recs
        if len(vcf_recs) == 0:
            stats["Pre-filter"]["No records after filter, discarded"] += 1
            write_to_log(log_hanlde, ss_rec_raw, "None; No records after filter; Discarded")
            continue

        # Discard ss_rec if there are multiple records
        if len(vcf_recs) > 1:
            stats["Pre-filter"]["Multiple records after filter, discarded"] += 1
            write_to_log(log_hanlde, ss_rec_raw, "None; >1 record after filter; Discarded")
            continue

        # Given that there is now 1 record, use that
        vcf_rec = vcf_recs[0]

        # Discard ssrec if there are multiple matching alleles
        if vcf_rec.n_alts() > 1:
            stats["Pre-filter"][">1 matching ref-alt pair for record, discarded"] += 1
            write_to_log(log_hanlde, ss_rec_raw, "None; >1 matching ref-alt pair for record; Discarded")
            continue

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
                    if (af_to_maf(ss_rec.eaf) > args.maf_palin_threshold or
                        af_to_maf(vcf_alt_af) > args.maf_palin_threshold):
                        stats["Palindromic"]["MAFs > maf_palin_threshold, discarded"] += 1
                        write_to_log(log_hanlde, ss_rec_raw, "Palindromic; MAFs > maf_palin_threshold; Discarded")
                        continue
                else:
                    stats["Palindromic"]["EAF not in sumstat file, discarded"] += 1
                    write_to_log(log_hanlde, ss_rec_raw, "Palindromic; EAF not in sumstat file; Discarded")
                    continue
                # Flip if eafs are not concordant
                if not afs_concordant(ss_rec.eaf, vcf_alt_af):
                    ss_rec.flip_beta()
                    ss_rec.revcomp_alleles()
                    stats["Palindromic"]["MAFs not concordant, alleles flipped"] += 1
                    write_to_log(log_hanlde, ss_rec_raw, "Palindromic; MAFs not concordant; Flipped")
                else:
                    stats["Palindromic"]["MAFs concordant, alleles correct"] += 1
                    write_to_log(log_hanlde, ss_rec_raw, "Palindromic; MAFs concordant; None")
            else:
                stats["Palindromic"]["infer_palin or infer_strand are False, discarded"] += 1
                write_to_log(log_hanlde, ss_rec_raw, "Palindromic; infer_palin or infer_strand False; Discarded")
                continue

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
                    stats["Reverse strand"]["alleles flipped"] += 1
                    write_to_log(log_hanlde, ss_rec_raw, "Reverse strand; alleles flipped; Flipped")
                else:
                    stats["Reverse strand"]["alleles correct"] += 1
                    write_to_log(log_hanlde, ss_rec_raw, "Reverse strand; alleles correct; None")
            else:
                stats["Reverse strand"]["infer_strand is False, discarded"] += 1
                write_to_log(log_hanlde, ss_rec_raw, "Reverse strand; infer_strand False; Discarded")
                continue

        # Harmonise same strand alleles
        elif compatible_alleles_forward_strand(ss_rec.other_al,
                                               ss_rec.effect_al,
                                               vcf_rec.ref_al,
                                               vcf_alt):
            # Flip if the effect allele matches the vcf ref alleles
            if ss_rec.effect_al.str() == vcf_rec.ref_al.str():
                ss_rec.flip_beta()
                stats["Forward strand"]["alleles flipped"] += 1
                write_to_log(log_hanlde, ss_rec_raw, "Forward strand; alleles flipped; Flipped")
            else:
                stats["Forward strand"]["alleles correct"] += 1
                write_to_log(log_hanlde, ss_rec_raw, "Forward strand; alleles correct; None")

        else:
            # Should never reach this 'else' statement
            sys.exit("Error: Alleles were not palindromic, opposite strand, or same strand!")

        #
        # Write ssrec to output ------------------------------------------------
        #

        # Add harmonised other allele, effect allele, eaf, beta to output
        out_row = ss_rec.data
        out_row["hm_hgvs"] = vcf_rec.hgvs()[0]
        out_row["hm_rsid"] = vcf_rec.id
        out_row["hm_other_allele"] = ss_rec.other_al.str()
        out_row["hm_effect_allele"] = ss_rec.effect_al.str()
        out_row["hm_eaf"] = str(ss_rec.eaf) if ss_rec.eaf else "NA"
        out_row["hm_beta"] = str(ss_rec.beta)

        # Write header
        if not header_written:
            outline = args.out_sep.join(out_row.keys()) + "\n"
            out_hanlde.write(outline.encode("utf-8"))
            header_written = True

        # Write row
        outline = args.out_sep.join(out_row.values()) + "\n"
        out_hanlde.write(outline.encode("utf-8"))

    # Print stats and write to log
    stat_str = process_stats_dict(stats)
    print("\n" + stat_str + "\n")
    outline = "\n# " + stat_str.replace("\n", "\n# ")
    log_hanlde.write(outline.encode("utf-8"))

    # Close handles
    out_hanlde.close()
    log_hanlde.close()

    print("Done!")

    return 0

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

    def hgvs(self):
        """ Represents variants in simplified HGVS format:
                chrom_pos_ref_alt
            Makes a list containing each alt allele
        Returns:
            list of strings
        """
        hgvs_list = []
        for alt in self.alt_als:
            hgvs = "{0}_{1}_{2}_{3}".format(self.chrom,
                                            self.pos,
                                            self.ref_al.str(),
                                            alt.str())
            hgvs_list.append(hgvs)
        return hgvs_list

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
        permissible = set([str(x) for x in list(range(1, 23)) + ["X", "Y", "MT"]])
        assert set([self.chrom]).issubset(permissible)
        # Assert that other and effect alleles are different
        assert self.other_al.str() != self.effect_al.str()

    def tolist(self):
        """ Returns identifying info (rsid, chrom, pos, alleles) as a list
        Returns:
            list
        """
        return [self.rsid, self.chrom, self.pos, self.other_al, self.effect_al]

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

def parse_info_field(field):
    """ Parses a vcf info field in to a dictionary
    Args:
        field (str): info field from vcf
    Returns:
        Dict: {name:[val1, val2]}
    """
    d = {}
    for entry in field.split(";"):
        try:
            key, value = entry.split("=")
            d[key] = value.split(",")
        except ValueError:
            d[entry] = []
    return d

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
        ss_record = SumStatRecord(row[args.rsid_col],
                                  row[args.chrom_col],
                                  row[args.pos_col],
                                  row[args.otherAl_col],
                                  row[args.effAl_col],
                                  row[args.beta_col],
                                  row.get(args.eaf_col, None),
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
    with open_gzip(inf, "rb") as in_handle:
        header = in_handle.readline().decode("utf-8").rstrip().split(sep)
        for line in in_handle:
            values = line.decode("utf-8").rstrip().split(sep)
            assert(len(values) == len(header))
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

def parse_args():
    """ Parse command line args using argparse.
    """
    parser = argparse.ArgumentParser(description="Summary statistc harmoniser")

    # Files
    parser.add_argument('--sumstats', metavar="<file>",
                        help=('GWAS summary statistics file'), type=str,
                        required=True)
    parser.add_argument('--vcf', metavar="<file>",
                        help=('Reference VCF file. Use # as chromosome wildcard.'), type=str, required=True)
    parser.add_argument('--out', metavar="<file>",
                        help=("Harmonised output file. Use 'gz' extension to gzip."), type=str,
                        required=True)
    parser.add_argument('--log', metavar="<file>",
                        help=("Log file. Use 'gz' extension to gzip."), type=str, required=True)
    # Columns
    parser.add_argument('--rsid_col', metavar="<str>",
                        help=('Rsid column'), type=str, required=True)
    parser.add_argument('--chrom_col', metavar="<str>",
                        help=('Chromosome column'), type=str, required=True)
    parser.add_argument('--pos_col', metavar="<str>",
                        help=('Position column'), type=str, required=True)
    parser.add_argument('--effAl_col', metavar="<str>",
                        help=('Effect allele column'), type=str, required=True)
    parser.add_argument('--otherAl_col', metavar="<str>",
                        help=('Other allele column'), type=str, required=True)
    parser.add_argument('--beta_col', metavar="<str>",
                        help=('beta column'), type=str, required=True)
    parser.add_argument('--eaf_col', metavar="<str>",
                        help=('EAF column'), type=str)
    # Other args
    parser.add_argument('--only_chrom', metavar="<str>",
                        help=('Only process provided chromosome.'), type=str)
    parser.add_argument('--maf_palin_threshold', metavar="<float>",
                        help=('Max MAF that will be used to infer palindrome strand (default: 0.42)'),
                        type=float, default=0.42)
    parser.add_argument('--af_vcf_min', metavar="<float>",
                        help=('Min freq of alt allele to be used (default: 0.001)'),
                        type=float, default=0.001)
    parser.add_argument('--af_vcf_field', metavar="<str>",
                        help=('VCF info field containing allele freq (default: AF_NFE)'),
                        type=str, default="AF_NFE")
    parser.add_argument('--in_sep', metavar="<str>",
                        help=('Input file column separator (default: tab)'),
                        type=str, default="\t")
    parser.add_argument('--out_sep', metavar="<str>",
                        help=('Output file column separator (default: tab)'),
                        type=str, default="\t")
    # Harmonisation options
    parser.add_argument('--infer_strand', metavar="<bool>",
                        help=('Infer and flip reverse strand alleles? [True|False] (default: True)'),
                        type=str2bool, default=True)
    parser.add_argument('--infer_palin', metavar="<bool>",
                        help=('Infer and flip palindromic alleles? [True|False] (default: True)'),
                        type=str2bool, default=True)
    return parser.parse_args()

def str2bool(v):
    """ Parses argpare boolean input
    """
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def process_stats_dict(stats):
    """ Process stats dictionary.
    Args:
        stats (dict of stats)
    Returns:
        str of
    """
    # Get number of records for each sub category
    sums = {}
    for key in stats:
        sums[key] = sum(stats[key].values())
    sums["Total"] = sum(sums.values())

    # Add record stats
    rows = ["Total records processed: {0}".format(sums["Total"])]
    for key in ["Pre-filter", "Palindromic", "Reverse strand", "Forward strand"]:
        perc = 100 * float(sums[key])/sums["Total"] if sums["Total"] != 0 else 0
        rows.append("{0} records: {1} ({2:.1f}%)".format(key, sums[key], perc))
        for key2 in stats[key]:
            perc = 100 * float(stats[key][key2])/sums["Total"] if sums["Total"] else 0
            rows.append("  {0}: {1} ({2:.1f}%)".format(key2, stats[key][key2], perc))

    return "\n".join(rows)

def initiate_stats():
    """ Returns dict of dicts which is used to store statistics
    """
    return {"Pre-filter":OrderedDict([
                           ("No record in VCF, discarded", 0),
                           ("No records after filter, discarded", 0),
                           ("Multiple records after filter, discarded", 0),
                           (">1 matching ref-alt pair for record, discarded", 0)
                           ]),
            "Palindromic":OrderedDict([
                           ("infer_palin or infer_strand are False, discarded", 0),
                           ("MAFs > maf_palin_threshold, discarded", 0),
                           ("EAF not in sumstat file, discarded", 0),
                           ("MAFs not concordant, alleles flipped", 0),
                           ("MAFs concordant, alleles correct", 0)
                           ]),
            "Reverse strand":OrderedDict([
                           ("alleles flipped", 0),
                           ("alleles correct", 0),
                           ("infer_strand is False, discarded", 0)
                           ]),
            "Forward strand":OrderedDict([
                           ("alleles flipped", 0),
                           ("alleles correct", 0)
                           ])
            }

def write_to_log(handle, ssrec, message):
    """ Write ssrecord to the log file with associated message.
    Args:
        handle (file): handle to log file
        ssrec (SumStatRecord): summary stat record
        message (str): message to log
    Returns: 0
    """
    outline = "\t".join([str(x) for x in ssrec.tolist()] + [message]) + "\n"
    handle.write(outline.encode("utf-8"))
    return 0

if __name__ == '__main__':

    main()
