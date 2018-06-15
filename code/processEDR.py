#
# Import necessary libraries
#
from datetime import datetime
from scipy import signal
import time
import pandas as pd
import numpy as np
import struct
import bitstring
from readLBL import parseLBLFile
from readAux import parseAuxFile
from readAnc import parseAncilliary
from readEDR import readEDRrecord, decompressSciData, detSCConf
from ProcessingTools import makeWindow
from SARProcessing import rangeCompression
from readChirp import detChirpFiles
from plottingFunctions import plotEDR_new, bytescl, rdr2san
import matplotlib.pyplot as plt
import glob, os, sys

def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')


def main(runName, auxname, lblname, edrname, presum_proc=None, beta=5, fil_type='Match', verb=True):
  """
    This python script simply reads in EDR data and returns the chirp
    compressed science record, which should be coplex voltage
  """
  logFile = '../runs/' + str(runName) + '.log'
  _log = open(logFile, 'w')
  writeLog(_log, 'Reading label file...')
  lblDic = parseLBLFile(lblname)
  writeLog(_log, 'Finished reading label file...')
  InstrPresum = lblDic['INSTR_MODE_ID']['Presum']       # Onboard presums
  if presum_proc is None:
    presum_proc = lblDic['INSTR_MODE_ID']['Presum']
  presum_fac = int(presum_proc / InstrPresum)
  if presum_fac < 1:
    presum_fac = 1
    print('WARNING: Processing presum less than onboard presumming....keeping native presum value') 
  instrMode = lblDic['INSTR_MODE_ID']['Mode']
  BitsPerSample = lblDic['INSTR_MODE_ID']['BitsPerSample']
  nrec = lblDic['FILE_RECORDS']                         # Number of records
  writeLog(_log, 'Reading Auxilliary file...')
  AuxDF = parseAuxFile(auxname, df=True)
  writeLog(_log, 'Finished reading Auxilliary file...')
  #
  # Determine to Bits per Sample
  #
  if BitsPerSample == 4:
    recLen = 1986
  elif BitsPerSample == 6:
    recLen = 2886
  elif BitsPerSample == 8:
    recLen = 3786
  if verb:
    writeLog(_log, 'Number of Records:\t{}'.format(nrec))
    writeLog(_log, '----- Presum Information -----')
    writeLog(_log, 'InstrPresum:\t{}'.format(InstrPresum))
    writeLog(_log, 'Processing Presum:\t{}'.format(presum_proc))
    if presum_fac == 1:
      writeLog(_log, 'No additional presumming performed')
    else:
      writeLog(_log, 'Presum Factor:\t{}'.format(presum_fac))
      writeLog(_log, '----- End Presum Information -----')
    writeLog(_log, 'Instrument Mode:\t{}'.format(instrMode))
    writeLog(_log, 'Bits Per Sample:\t{}'.format(BitsPerSample))
    writeLog(_log, 'Record Length:\t{}'.format(recLen))
  #
  #
  # Construct window
  #
  window, win_str = makeWindow(beta, length=2048)
  writeLog(_log, win_str)
  ######################################################################
  #
  # Begin Processing
  #
  ######################################################################
  EDRData = np.zeros([2048, int(np.ceil(nrec/presum_fac))], complex)
  presum_rec = np.zeros([2048, presum_fac], complex)
  writeLog(_log, 'Opening EDR File:\t{}'.format(edrname))
  _file1 = open(edrname, 'rb')
  if verb:
    writeLog(_log, 'Begin decompression at:\t{}'.format(datetime.now()))
  for _i in range(0, nrec, presum_fac):
    for _k in range(0, presum_fac):
      #
      # Get the time since the start of the run for this pulse group
      #
      tcen = AuxDF['ELAPSED_TIME'][_i]
      #
      # Read in single record
      #
      sci, ancil = readEDRrecord(_file1, _i, recLen, BitsPerSample)
      #
      # Parse Ancilliary Datq
      #
      ancil = parseAncilliary(ancil)
      #
      # Decompress science data
      #
      sci = decompressSciData(sci, lblDic['COMPRESSION'], InstrPresum, BitsPerSample, ancil['SDI_BIT_FIELD'])
      #
      # Determine calibrated chirp
      #
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i])
      #
      # Perform chirp compression
      #
      presum_rec[:,_k] = rangeCompression(sci, calChirp, window, fil_type="Match", diag=False)
      #
      # Perform on-ground calibration
      #  This step will have to wait. Apparently the angles given in the Auxilliary file are
      #  incorrect for some dates and change after a certain date. Fabrizio is checking on this
      #
#      presum_rec[:, _k] = calibrateData(presum_rec[:, _k])
      if _k == 0:
        mx0 = np.where(presum_rec[:,_k] == np.max(presum_rec[:,_k]))[0]
      else:
        mx = np.where(presum_rec[:,_k] == np.max(presum_rec[:,_k]))[0]
        sft = mx0 - mx
        presum_rec[:,_k] = np. roll(presum_rec[:, _k], sft)
    EDRData[:,int(_i/presum_fac)] = np.sum(presum_rec, axis=1)
    #temp = np.power(np.abs(EDRData[:,int(_i/presum_fac)]),2)/ np.max(np.power(np.abs(EDRData[:,int(_i/presum_fac)]),2))
  '''
    temp = np.abs(np.real(EDRData[:,int(_i/presum_fac)]))
    temp = temp / np.max(temp)
    idx = np.where(temp == np.max(temp))[0][0]
    plt.plot(temp)
    plt.xlim(idx, idx+100)
    plt.show()
    if _i == 5*presum_fac:
      sys.exit()
  '''
  if verb:
    writeLog(_log, 'Decompession finished at:\t{}'.format(datetime.now()))
  fname = '../runs/' + str(runName) + '.npy'
  np.save(fname, EDRData)
  #plotEDR(EDRData, fname=runName, rel=True)
  plotEDR_new(EDRData, fname=runName)
  return
  
if __name__ == '__main__':
  runName = 'td7_beta0_ps4'
  verb = True
  win_type = 14                                         # 0 (uniform), 2 (bartlett), 3 (Hann), 4 (Hamming), 5 (Blackman), 6 (Kaiser)
  beta = 0
  td = 7                                               # Which test set
  fil_type = 'Match'                                    # Chirp compression method
  presum_proc = 4
  #
  # BEGIN INPUT DATA
  #
  if td == 0:
    auxname = '../data/e_0168901_001_ss19_700_a_a.dat'
    lblname = '../data/e_0168901_001_ss19_700_a.lbl'
    edrname = '../data/e_0168901_001_ss19_700_a_s.dat'
  elif td == 1:
    auxname = '../data/e_1601301_001_ss19_700_a_a.dat'
    lblname = '../data/e_1601301_001_ss19_700_a.lbl'
    edrname = '../data/e_1601301_001_ss19_700_a_s.dat'
  elif td == 2:
    auxname = '../data/e_1408801_001_ss11_700_a_a.dat'
    lblname = '../data/e_1408801_001_ss11_700_a.lbl'
    edrname = '../data/e_1408801_001_ss11_700_a_s.dat'
  elif td == 3:
    auxname = '../data/e_0169401_001_ss18_700_a_a.dat'
    lblname = '../data/e_0169401_001_ss18_700_a.lbl'
    edrname = '../data/e_0169401_001_ss18_700_a_s.dat'
  elif td == 4: #VERY FLAT REGION
    auxname = '../data/e_0323403_001_ss19_700_a_a.dat'
    lblname = '../data/e_0323403_001_ss19_700_a.lbl'
    edrname = '../data/e_0323403_001_ss19_700_a_s.dat'
  elif td == 5: #ANOTHER SHORT SEMI FLAT REGION
    auxname = '../data/e_2777101_001_ss19_700_a_a.dat'
    lblname = '../data/e_2777101_001_ss19_700_a.lbl'
    edrname = '../data/e_2777101_001_ss19_700_a_s.dat'
  elif td == 6: #Bruce's test region
    auxname = '../data/e_0879801_001_ss19_700_a_a.dat'
    lblname = '../data/e_0879801_001_ss19_700_a.lbl'
    edrname = '../data/e_0879801_001_ss19_700_a_s.dat'
  elif td == 7: #NPLD Test; many reflections
    auxname = '../data/e_0577001_001_ss19_700_a_a.dat'
    lblname = '../data/e_0577001_001_ss19_700_a.lbl'
    edrname = '../data/e_0577001_001_ss19_700_a_s.dat'
  elif td == 8: #Rolled Observation
    auxname = '../data/e_4866701_001_ss19_700_a_a.dat'
    lblname = '../data/e_4866701_001_ss19_700_a.lbl'
    edrname = '../data/e_4866701_001_ss19_700_a_s.dat'
  elif td == 9: #Rolled Observation
    auxname = '../data/e_1733901_001_ss19_700_a_a.dat'
    lblname = '../data/e_1733901_001_ss19_700_a.lbl'
    edrname = '../data/e_1733901_001_ss19_700_a_s.dat'
  elif td == 10: # looooooong obs.
    auxname = '../data/e_0246001_001_ss11_700_a_a.dat'
    lblname = '../data/e_0246001_001_ss11_700_a.lbl'
    edrname = '../data/e_0246001_001_ss11_700_a_s.dat'
  elif td == 11: # looooooong obs.
    auxname = '../data/e_5050702_001_ss19_700_a_a.dat'
    lblname = '../data/e_5050702_001_ss19_700_a.lbl'
    edrname = '../data/e_5050702_001_ss19_700_a_s.dat'
     
  main(runName, auxname, lblname, edrname, presum_proc=presum_proc, beta=beta, fil_type=fil_type, verb=verb) 
