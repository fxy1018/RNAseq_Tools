'''
Created on Jul 20, 2017

@author: XFan
'''

#this is a command tool to run quantifier rna-seq tools
import os
import logging
from FileDealer import FileObj
import rpy2


class Quantifier(object):
    def __init__(self):
        pass
    
    def countRead(self, bamFile, gtfFile, outFile): 
        #use version 0.9.0, bioconda only have version 0.7
        cmd = "python -m HTSeq.scripts.count -f bam --additional-attr=gene_name --max-reads-in-buffer=30000000" + \
            " " + bamFile.path + " " + gtfFile.path + " >> " + outFile.path
        
        try:
#             print(cmd)
            os.system(cmd)
        except:
            mss = bamFile.fileName + " failed to executive"
            logging.warning(mss)
 
