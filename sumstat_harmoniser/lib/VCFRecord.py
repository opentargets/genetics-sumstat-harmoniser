from lib.Seq import Seq

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
        try:
            self.info = parse_info_field(row[7])
        except IndexError:
            self.info = {}

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
