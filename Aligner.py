'''
Created on Jul 18, 2017

@author: XFan
'''


'''

This is for running aligner on cmd on linux

'''
from FileDealer import FileObj
import os
import logging


class Aligners(object):
    def __init__(self):
        pass
    
    def STAR(self, index_dir, read, gtf, out_dir):
        if not os.path.exists(out_dir.dir):
            os.system("sudo mkdir " + out_dir.dir )
             
        cmd = "sudo /home/xfan/anaconda3/envs/aligners/bin/STAR" + \
             " --runThreadN 2"  + \
             " --genomeDir " + index_dir + \
             " --readFilesIn " + read.path + \
             " --readFilesCommand sudo gunzip -c" + \
             " --sjdbGTFfile " + gtf.path + \
             " --sjdbOverhang 75" + \
             " --outFileNamePrefix " + out_dir.path + \
             " --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted"
        
        zip_cmd = "sudo gzip " + out_dir.path + "Aligned.out.bam"
        print(cmd)
        print(zip_cmd)
        try:
            os.system(cmd)
            os.system(zip_cmd)
        except:
            mss = self.file + " failed to executive"
            logging.warning(mss)
          
        
    def tophat(self):
        pass
     
    def novoalign(self):
        pass
    
def getAlignReport(dirs):
    pass
     