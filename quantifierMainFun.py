'''
Created on Jul 24, 2017

@author: XFan

'''
from Quantifier import Quantifier
from FileDealer import FileObj
import multiprocessing


def fun(b, gtfFile):
    q = Quantifier()
    name = b.fileName.strip("Aligned.out.bam")
    name += "_count.txt"
    outFile = FileObj(directory="/home/xfan/htseqTest/", name=name)
    q.countRead(b, gtfFile, outFile)

if __name__ == "__main__":
    
    
    fastqs = ["3.fastq", "4.fastq", "5.fastq", "6.fastq", "7.fastq", "8.fastq", "9.fastq", "10.fastq", "11.fastq","12.fastq"]
#     fastqs = ["2.fastq"]
    bamroot = "/data/cofactor_genomics/IWP0005JJ/alignment/genome/"
    bamFiles = []
    
    gtfFile = FileObj(directory="/data/cofactor_genomics/reference/annotation/", name = "Homo_sapiens.GRCh38.89.gtf")
    
    for f in fastqs:
        dir = bamroot + f + "_aligned" 
        fileName = f + "Aligned.out.bam"
        fileObj = FileObj(dir, fileName)
        bamFiles.append(fileObj)
    
    jobs = []
    for b in bamFiles:
        p = multiprocessing.Process(target = fun, args = (b, gtfFile))
        jobs.append(p)
        p.start()
        
        

