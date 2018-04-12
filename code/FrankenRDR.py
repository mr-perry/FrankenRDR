#
# Import Libraries
#
from readLBL import loadLBL, parseLBL
from readAuxFile import readAuxFile, readRow, detCalibChirp, loadCalChirp
#from readEDR import parseEDRname, loadData
import matplotlib.pyplot as plt
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
  labelFile = [s for s in files if ".lbl" in s][0]
  scienceFile = [s for s in files if 's.dat' in s][0]
  auxFile = [s for s in files if 'a.dat' in s][0]
  #
  # Parse label files
  #
  lblFile = loadLBL(labelFile)
  lblDict = parseLBL(lblFile)
  print(lblDict)
  sys.exit()
  #
  # Load auxillary data
  # 
  _auxFile = readAuxFile(auxFile) 
  #
  # This little blurb was for testing purposes
  #
  '''
  rowNum = range(6,7)
  if _auxFile != []:
    for _row in rowNum:
      AuxData = readRow(auxFile, _row)
      print(AuxData['TX_TEMP'], AuxData['RX_TEMP'])
      calChirpFile = detCalibChirp(AuxData['TX_TEMP'], AuxData['RX_TEMP'])
      calChirp, calChirp_reorg = loadCalChirp(calChirpFile)
      plt.subplot(2,1,1)
      plt.plot(calChirp)
      plt.subplot(2,1,2)
      plt.plot(calChirp_reorg)
      plt.show() 
  sys.exit()
  '''
  #
  # Parse EDR file name
  #
  #print(scienceFile)
  #EDRData = parseEDRname(scienceFile)
  #print(EDRData)
  return

if __name__ == '__main__':
  main()
