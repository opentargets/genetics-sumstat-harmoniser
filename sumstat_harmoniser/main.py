#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#
# Version 0.2
#
# Harmonise GWAS summary statistics against a reference VCF

import os
import sys
import gzip
from copy import deepcopy
from subprocess import Popen, PIPE
from collections import OrderedDict, Counter
from lib.VCFRecord import VCFRecord
from lib.SumStatRecord import SumStatRecord
from lib.Seq import Seq
from lib.functions import *

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


if __name__ == '__main__':

    main()
