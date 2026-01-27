####bin/expAnalysis.py
import subprocess
import sys
import os

def expAnalysis(pepQuant,transcriptQuant,cor,outpath,manifest=None):
  pepQuant=os.path.abspath(pepQuant)
  transcriptQuant=os.path.abspath(transcriptQuant)
  if not outpath:
    outpath='/'
    outpath=outpath.join(pepQuant.split('/')[:-1])+'/expressionAnalysis'
    if not os.path.exists(outpath):
      os.makedirs(outpath)
  if not manifest:
    manifest=''
  scriptPath=sys.path[0]+'/expAnalysis/'
  process_cmd='Rscript '+scriptPath+'expAnalysis.R -p '+pepQuant+' -t '+transcriptQuant+' -c '+cor+' -m '+manifest+' -o '+outpath
  subprocess.run(process_cmd,shell=True,check=True)

if __name__ == '__main__':
    pepQuant = "/home/yvzeng/Project/software/proteomics/test/quant/quant.tsv"
    transcriptQuant = "/home/yvzeng/Project/software/proteomics/test/quant/transcriptQuantTest.xls"
