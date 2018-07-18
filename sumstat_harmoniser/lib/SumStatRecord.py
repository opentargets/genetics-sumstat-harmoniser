from lib.Seq import Seq

class SumStatRecord:
    """ Class to hold a summary statistic record.
    """
    def __init__(self, chrom, pos, other_al, effect_al, beta, oddsr, eaf,
                 data):

        # Set raw info
        self.chrom = str(chrom)
        self.pos = int(pos)
        self.other_al = Seq(other_al)
        self.effect_al = Seq(effect_al)
        self.data = data
        self.beta = float(beta) if beta else None
        self.oddsr = float(oddsr) if oddsr else None

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

        # Set harmonised values
        self.hm_rsid = None
        self.hm_chrom = None
        self.hm_pos = None
        self.hm_other_al = None
        self.hm_effect_al = None
        self.is_harmonised = False
        self.hm_code = None

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
