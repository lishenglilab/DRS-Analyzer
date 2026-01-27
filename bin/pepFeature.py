####bin/pepFeature.py
import subprocess
import sys
import os

def pepFeature(quantFile,pepRef,transcriptRef,outpath):
  quantFile=os.path.abspath(quantFile)
  pepRef=os.path.abspath(pepRef)
  pepRef=os.path.abspath(pepRef)
  if not outpath:
    outpath='/'
    outpath=outpath.join(quantFile.split('/')[:-1])+'/pepFeature'
    if not os.path.exists(outpath):
      os.makedirs(outpath)
  scriptPath=sys.path[0]+'/pepFeature/'
  process_cmd='Rscript '+scriptPath+'pepFeature.R -q '+quantFile+' -p '+pepRef+' -t '+transcriptRef+' -o '+outpath
  subprocess.run(process_cmd,shell=True,check=True)
