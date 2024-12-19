class Seq:
    @staticmethod
    def reverse_complement(seq):
        """
        param DNA sequence: eg ATAG
        return reverse_complement sequence : CTAT
        """
        seq = seq[::-1].upper()
        return seq.replace('A', 't').\
                    replace('C', 'g').\
                    replace('T', 'a').\
                    replace('G', 'c').upper()
    def GC_content(seq):
        """
        param DNA sequence: eg ATAG
        return GC content : 0.5
        method : (G+C)/(G+C+A+T)
        """
        seq = seq.upper()
        return (seq.count('G') + seq.count('C')) / len(seq)
    def GC_skew(seq):
        """
        param DNA sequence: eg ATAG
        return GC skew : 0.5
        method : (G-C)/(G+C)
        """
        return (seq.count('G') - seq.count('C')) /  (seq.count('G') + seq.count('C'))


acid_info = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Ter': '*'
}
#@function_timer
def get_oneletter_hgvsp(hgvsp):
    '''获取单字母的hgvsp结果
    '''
    for acid in acid_info:
        if acid in hgvsp:
            hgvsp = hgvsp.replace(acid, acid_info[acid])
    return hgvsp
