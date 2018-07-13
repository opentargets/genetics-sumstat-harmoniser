import argparse

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
    # parser.add_argument('--infer_strand', metavar="<bool>",
    #                     help=('Infer and flip reverse strand alleles? [True|False] (default: True)'),
    #                     type=str2bool, default=True)
    parser.add_argument('--palin_mode', metavar="[infer|forward|reverse|drop]",
                        help=('Mode to use for palindromic variants: '
                              '(a) infer strand from effect-allele freq, '
                              '(b) assume forward strand, '
                              '(c) assume reverse strand, '
                              '(d) drop all palindromic variants. '
                              '(default: infer)'),
                        choices=['infer', 'forward', 'reverse', 'drop'],
                        type=str2bool, default='infer')
    return parser.parse_args()

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
