#
# Import necessary libraries
#
from datetime import datetime
from scipy.signal import welch
import time
import pandas as pd
import numpy as np
import struct
import bitstring
from readLBL import parseLBLFile
from readAux import parseAuxFile
from readAnc import parseAncilliary
from readEDR import readEDRrecord, decompressSciData, detSCConf
from ProcessingTools import makeWindow, calcSNR
from SARProcessing import rangeCompression
from readChirp import detChirpFiles
from plottingFunctions import plotEDR, bytescl, rdr2san, makePSD, plotFirstReturn 
import matplotlib.pyplot as plt
import glob, os, sys

def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')


def main(runName, auxname, lblname, edrname, chirp='ref', presum_proc=None, beta=5, fil_type='Match', verb=True, diag=False):
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
  if chirp == 'ideal' or chirp == 'UPB':
    window, win_str = makeWindow(beta, length=3600)
  else:
    window, win_str = makeWindow(beta, length=2048)
  writeLog(_log, win_str)
  ######################################################################
  #
  # Begin Processing
  #
  ######################################################################
  if chirp == 'ideal' or chirp == 'UPB':
    if chirp == 'ideal':
      print('Using ideal chirp')
    else:
      print('Using cal_filter chirp from UPB')
    EDRData = np.zeros([3600, int(np.ceil(nrec/presum_fac))], complex)
    presum_rec = np.zeros([3600, presum_fac], complex)
  else:
    EDRData = np.zeros([4096, int(np.ceil(nrec/presum_fac))], complex)
    presum_rec = np.zeros([4096, presum_fac], complex)
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
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i], chirp=chirp)
      #
      # Perform chirp compression
      #
      presum_rec[:,_k] = rangeCompression(sci, calChirp, window, chirp=chirp, fil_type="Match", diag=diag)
      #
      # Perform on-ground calibration
      #  This step will have to wait. Apparently the angles given in the Auxilliary file are
      #  incorrect for some dates and change after a certain date. Fabrizio is checking on this
      #
#      presum_rec[:, _k] = calibrateData(presum_rec[:, _k])
    EDRData[:,int(_i/presum_fac)] = np.sum(presum_rec, axis=1)
    plotFirstReturn(EDRData[:,int(_i/presum_fac)], type='Amp', sidelobe=False, title='First Return', fname='VibraSeis', dpi=500)
    sys.exit()
  #
  # Calculate Signal-To-Noise (This isn't very useful for determining the best method for reducing sidelobe
  #
  SNR = calcSNR(EDRData)
  if verb:
    writeLog(_log, 'Decompession finished at:\t{}'.format(datetime.now()))
    writeLog(_log, 'Signal-To-Noise Ratio:\t{}'.format(SNR))
  fname = '../runs/' + str(runName) + '.npy'
  np.save(fname, EDRData)
  plotEDR(EDRData, fname=runName, ptype='Amp', thres=0, rel=True)
  plotEDR(EDRData, fname=runName, ptype='Pow', thres=0, rel=True)
  plotEDR(EDRData, fname=runName, ptype='dB', thres=-15, rel=True)
  return
  
if __name__ == '__main__':
  runName = 'test'
  verb = True
  diag = True
  chirp = 'ref'
  beta = 0						# Kaiser window beta value; 0 -rectangular; 5 similar to Hamming; 6 similar to Hann, 8.6 Similar to blackman
  td = 6                                              # Which test set
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
  elif td == 8: 
    auxname = '../data/e_3940901_001_ss19_700_a_a.dat'
    lblname = '../data/e_3940901_001_ss19_700_a.lbl'
    edrname = '../data/e_3940901_001_ss19_700_a_s.dat'
     
  main(runName, auxname, lblname, edrname, chirp=chirp, presum_proc=presum_proc, beta=beta, fil_type=fil_type, verb=verb, diag=diag) 
