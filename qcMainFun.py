'''
Created on Jul 14, 2017

@author: XFan
'''
from FileDealer import FileObj
from FileDealer import getSeqFiles
from FileDealer import makeFastDirs
from QualityController import QualityController
import os
from QualityController import fastqcReport

 
if __name__ == "__main__" :
    #build the folders for all samples
#     expNum = ["IWP0002JJ","IWP0003JJ","IWP0004JJ","IWP0005JJ"]
#     data_root = '/data/cofactor_genomics/'q
#     fastqc_root = '/home/xfan/qcTest'
#     
#     for e in expNum:
#         seqFiles = getSeqFiles(data_root, e)
#     
#         for i in seqFiles:
#             makeFastDirs(fastqc_root, e, i)

    #view quality of original seq data 
#     expNum = "IWP0005JJ"
#     data_root = '/data/cofactor_genomics/'
#     fastqc_root = '/home/xfan/qcTest'
#     seqFiles = getSeqFiles(data_root, expNum)
#     
#     for i in seqFiles:
#         data_path = os.path.join(data_root, expNum, i)
#         qcResult_path = os.path.join(fastqc_root, expNum, i, "raw")
#         qcControl = QualityController()
#         qcControl.FastQC(data_path, qcResult_path)
#     

    #view quality of original seq data 
#     expNum = "IWP0005JJ"
#     data_root = '/data/cofactor_genomics/'
#     fastqc_root = '/home/xfan/qcTest'
#     seqFiles = getSeqFiles(data_root, expNum)
#     qcReport(fastqc_root, expNum, seqFiles)   


#to filter out rRNA from raw reads
    qc = QualityController()
    ref_root = "/data/cofactor_genomics/reference/rRNA"
    refs = ["rfam_5s_id98.fasta", "rfam_5.8s_id98.fasta", "silva_euk_28s_id98.fasta", "silva_euk_18s_id95.fasta"]
    refObjs = []
    for r in refs:
        file = FileObj(ref_root, r)
        refObjs.append(file)
    
    index_root = "/data/cofactor_genomics/reference/rRNA/sortmeRNA_index"
    indexs = ["rfam_5s_id98_db", "rfam_5.8s_id98_db", "silva_euk_28s_id98_db", "silva_euk_18s_id95_db"]
    indexObjs = []
    for i in indexs:
        file = FileObj(index_root, i)
        indexObjs.append(file)
    

    expNum = "IWP0005JJ"
    data_root = '/data/cofactor_genomics/'
    seq_files = getSeqFiles(data_root, expNum)
    for f in seq_files:
        readObj = FileObj("/data/cofactor_genomics/IWP0005JJ", f)
        qc.sortmerna(refObjs, indexObjs, readObj)

        
    
    
    
    
    
    