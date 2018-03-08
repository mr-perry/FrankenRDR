#
# Import Libraries
#
from readAuxFile import readAuxFile, readRow, detCalibChirp, loadCalChirp
#from readEDR import parseEDRname, loadData
import glob, os

def main():
  OrbitNum = 1689
  dataPath = '../data/'
  OrbitNum = str(OrbitNum)
  #
  # Check length and pad the front with a zero if necessary
  #
  if len(OrbitNum) == 4:
    OrbitNum = '0' + OrbitNum
  #
  # Find data associated with orbit number
  #
  fileSearch = dataPath + '?_'+OrbitNum+'??_*'
  files = glob.glob(fileSearch)
  print(files)
  labelFile = [s for s in files if ".lbl" in s][0]
  scienceFile = [s for s in files if 's.dat' in s][0]
  auxFile = [s for s in files if 'a.dat' in s][0]
  #
  # Parse label files
  #
#  labels = parseLabelFile()
  #
  # Load science data
  # 
  print(scienceFile)
  print(auxFile)
  return

if __name__ == '__main__':
  main()
  
   
