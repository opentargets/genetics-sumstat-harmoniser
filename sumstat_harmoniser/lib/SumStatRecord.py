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
