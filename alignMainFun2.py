'''
Created on Jul 21, 2017

@author: XFan
'''


from FileDealer import FileObj
from FileDealer import getSeqFiles
from FileDealer import makeFastDirs
from FileDealer import sam2bam
from Aligner import Aligners
import os

if __name__ == "__main__":
    aligner = Aligners()
    index_dir = "/data/cofactor_genomics/reference/hg38/STAR_index"
    gtf = FileObj("/data/cofactor_genomics/reference/annotation", "Homo_sapiens.GRCh38.89.gtf")
#     base_files = ["5.fastq", "6.fastq", "7.fastq", "8.fastq", "9.fastq"]
    base_files = ["4.fastq"]
  
    for file in base_files:
        read_name = file.strip(".fastq")
        read_name += "_non_rRNA.fastq.gz"
        read = FileObj("/data/cofactor_genomics/IWP0005JJ/rRNA_filtered", read_name)
        dir_name = file + "_aligned"
        out_dir = FileObj("/data/cofactor_genomics/IWP0005JJ/alignment/genome/" + dir_name + "/", file)
        aligner.STAR(index_dir, read, gtf, out_dir)

        
        
    

 