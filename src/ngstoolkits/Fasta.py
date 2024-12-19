import pysam
from .sequence import SeqFunction

class Fasta():
    """
    class for fasta file,
    with some useful function when dealing with fasta file
    """
    def __init__(self, fasta:str, age):
        self.reference = pysam.FastaFile(fasta)

    def gc_rate_special_region(self,chrom:str,start:int,end:int):
        """
        param DNA sequence: eg ATAG
        return GC content : 0.5
        method : (G+C)/(G+C+A+T)
        """
        seq = self.reference.fetch(chrom, start, end)
        seq = seq.upper()
        return SeqFunction.GC_content(seq)
