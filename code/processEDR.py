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
from readEDR import readEDRrecord, decompressSciData
from ProcessingTools import makeWindow
from SARProcessing import rangeCompression
from readChirp import detChirpFiles
from plottingFunctions import plotEDR, bytescl, rdr2san
import matplotlib.pyplot as plt
import glob, os, sys

def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')


def main(runName, auxname, lblname, edrname, win_type=4, fil_type='Match', compression=False, verb=True):
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
    writeLog(_log, 'InstrPresum:\t{}'.format(InstrPresum))
    writeLog(_log, 'Instrument Mode:\t{}'.format(instrMode))
    writeLog(_log, 'Bits Per Sample:\t{}'.format(BitsPerSample))
    writeLog(_log, 'Record Length:\t{}'.format(recLen))
  #
  #
  # Construct window
  #
  window, win_str = makeWindow(win_type)
  writeLog(_log, win_str)
  ######################################################################
  #
  # Begin Processing
  #
  ######################################################################
  EDRData = np.zeros([2048, nrec], complex)
  writeLog(_log, 'Opening EDR File:\t{}'.format(edrname))
  _file1 = open(edrname, 'rb')
  if verb:
    writeLog(_log, 'Begin decompression at:\t{}'.format(datetime.now()))
  for _i in range(0, nrec):
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
    sci = decompressSciData(sci, compression, InstrPresum, BitsPerSample)
    #
    # Determine calibrated chirp
    #
    calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i])
    #
    # Perform chirp compression
    #
    EDRData[:,_i] = rangeCompression(sci, calChirp, window, fil_type="Match", diag=False)
  if verb:
    writeLog(_log, 'Decompession finished at:\t{}'.format(datetime.now()))
    writeLog(_log, 'Converting to presum:\t{}'.format(presum_proc))
  fname = '../runs/' + str(runName) + '.npy'
  np.save(fname, EDRData)
  plotEDR(EDRData, fname=runName, rel=True)
  return
  
if __name__ == '__main__':
  runName = 'Test'
  verb = True
  win_type = 4
  td = 6						# Which test set
  fil_type = 'Match'					# Chirp compression method
  compression = False;					# Static compression
  presum_proc = 16					# Bring processing to presum-16
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
  main(runName, auxname, lblname, edrname, win_type=win_type, fil_type=fil_type, compression=compression, verb=verb) 
