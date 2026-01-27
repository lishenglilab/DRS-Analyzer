####bin/dataClean.py
import subprocess
import sys
import os

def dataCleaning(quantFile):
  quantFile=os.path.abspath(quantFile)
  scriptPath=sys.path[0]+'/dataCleaning/'
  process_cmd='Rscript '+scriptPath+'dataProcessing.R -q '+quantFile
  subprocess.run(process_cmd,shell=True,check=True)
