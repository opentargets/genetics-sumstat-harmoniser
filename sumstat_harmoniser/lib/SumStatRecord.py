from lib.Seq import Seq

class SumStatRecord:
    """ Class to hold a summary statistic record.
    """
    def __init__(self, chrom, pos, other_al, effect_al, beta, oddsr,
                 oddsr_lower, oddsr_upper, eaf, data):

        # Set raw info
        self.chrom = chrom
        self.pos = pos
        self.other_al = other_al
        self.effect_al = effect_al
        self.data = data
        self.beta = float(beta) if isFloat(beta) else None
        self.oddsr = float(oddsr) if isFloat(oddsr) else None
        self.oddsr_lower = float(oddsr_lower) if isFloat(oddsr_lower) else None
        self.oddsr_upper = float(oddsr_upper) if isFloat(oddsr_upper) else None

        # Effect allele frequency is not required if we assume +ve strand
        if isFloat(eaf):
            self.eaf = float(eaf)
            assert 0<= self.eaf <= 1
        else:
            self.eaf = None

        # Set harmonised values
        self.hm_rsid = None
        self.hm_chrom = None
        self.hm_pos = None
        self.hm_other_al = None
        self.hm_effect_al = None
        self.is_harmonised = False
        self.hm_code = None

    def validate_ssrec(self):
        ''' Ensures that chrom, pos, other_al, effect_al are of correct type
            Return code which will either be:
                - None if successful,
                - 14 if unsuccessful
        '''
        # Coerce types
        self.chrom = str(self.chrom)
        self.other_al = Seq(self.other_al)
        self.effect_al = Seq(self.effect_al)
        try:
            self.pos = int(self.pos)
        except (ValueError, TypeError) as e:
            return 14

        # Assert that other and effect alleles are different
        if self.other_al.str() == self.effect_al.str():
            return 14

        # Assert that unkown nucleotides don't exist
        if 'N' in self.other_al.str() or 'N' in self.effect_al.str():
            return 14

        return None

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
        if self.beta:
            self.beta = self.beta * -1
        # Flip OR
        if self.oddsr:
            self.oddsr = self.oddsr ** -1
        if self.oddsr_lower and self.oddsr_upper:
            unharmonised_lower = self.oddsr_lower
            unharmonised_upper = self.oddsr_upper
            self.oddsr_lower = unharmonised_upper ** -1
            self.oddsr_upper = unharmonised_lower ** -1
        # Switch alleles
        new_effect = self.other_al
        new_other = self.effect_al
        self.other_al = new_other
        self.effect_al = new_effect
        # Flip eaf
        if self.eaf:
            self.eaf = 1 - self.eaf

    def alleles(self):
        """
        Returns:
            Tuple of (other, effect) alleles
        """
        return (self.other_al, self.effect_al)

    def __repr__(self):
        return "\n".join(["Sum stat record",
                          "  chrom        : " + self.chrom,
                          "  pos          : " + str(self.pos),
                          "  other allele : " + str(self.other_al),
                          "  effect allele: " + str(self.effect_al),
                          "  beta         : " + str(self.beta),
                          "  odds ratio   : " + str(self.oddsr),
                          "  EAF          : " + str(self.eaf)
                          ])


def isFloat(value):
    if value is not None:
        try:
            float(value)
            return True
        except ValueError:
            return False
    else:
        return False
