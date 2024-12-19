import pysam

class BAM:
    def __init__(self,Bam:pysam.AlignmentFile):
        self.Bam = Bam
    def get_base_count_per_position(self,chr:str,pos:int):
        """
        Get the base count per position
        :param Bam: pysam.AlignmentFile
        :param chr: str
        :param pos: int
        :return: numbers(int) of each base in order [A C G T total] at the position
        """
        count=self.Bam.count_coverage(chr, pos-1,pos,quality_threshold=0)
        return [count[0][0],count[1][0],count[2][0],count[3][0],count[0][0]+count[1][0]+count[2][0]+count[3][0]]

    def get_region_depth(self,chr:str,start:int,end:int=0,flank:int=10):
        """
        Get the average depth at the position
        :param Bam: pysam.AlignmentFile
        :param chr: str
        :param pos: int
        :param window_size: int
        :return: average_depth, region_detail
        average_depth: average depth at the position
        region_detail: list of base count per position
        """
        if end==0:
            Warning("end is not set, use flank:" + flank)
            end=start+flank
            start=start-flank
        elif start>end:
            Warning("start is larger than end, change start and end")
            start,end=end,start
        region_site_num=0
        region_cover_base=0
        region_detail=[]
        for position in range(start,end):
            region_site_num+=1
            region_cover_base+=self.get_base_count_per_position(self.Bam,chr,position)[4]
            region_detail.append(self.get_base_count_per_position(self.Bam,chr,position))
        region_average_depth=region_cover_base/region_site_num
        return region_average_depth,region_detail
