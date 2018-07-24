class Seq:
    """ DNA nucleotide class
    Args:
        seq (str): dna sequence
    """
    def __init__(self, seq):
        self.seq = str(seq).upper()
        assert set(list(self.seq)).issubset(set(["A", "T", "G", "C", "D", "I", "N", "."]))

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
