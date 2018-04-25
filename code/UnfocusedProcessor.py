#
# Import Libraries
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
from readChirp import detChirpFiles
from makeWindow import makeWindow
import matplotlib.pyplot as plt
import glob, os, sys



def writeLog(log, string):
  print(string)
  if log is not None:
    log.write(string + '\n')

def main():
  """
    This is an attempt to recreate the Italian RDR processor in the style of the unfocused processor
    created by BC. The script will be based off reading data as it is set up in the PDS (i.e. *a_a.dat,
    *a_s.dat, and *.lbl files).
    
    User options (these are currently set in the code itself):
      verb: Print processing info (Default: True)
      win_type: Set the window type (Default: 1; cosine bell)
  """
  global nrec, instrMode, InstrPresum, BitsPerSample, AuxDF, logFile 
  verb = True
  win_type = 3
  minsundist = 206644545
  delay_res = 135.00e-06 / 3600.		# Delay resolution (sec) of each bin
  sharad_ipp = 1. / 700.28			# Time (sec) between start of each pulse
  c = 3.000e08
  ###################################################################################
  # 
  # TEST DATA
  #
  td = 4					# Select Test data set (0 - 4)
  fil_type = 'Match'					# Method of filtering the chirp (Match or Inverse
  fft_length = 1536
  compression = False;
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
    ncols = 448						# These need to be hardcoded until I can determine
  #
  # Define output files
  #
  file2 = '../runs/unfocgram.raw'			# 32-bit file
  file2a = '../runs/unfocpow.raw'			# 32-bit power
  file3 = '../runs/unfocgram.tif'			# RGB tiff
  file3a = '../runs/unfocpow.tif'			# Power tiff
  cplxfile = '../runs/complexgram.raw'			# Output complex file
  ###############################################################################################
  logFile = '../runs/log.txt'
  #
  # Open Log file
  #
  _log = open(logFile, 'w')
  ###############################################################################################
  writeLog(_log, 'Reading label file...')
  lblDic = parseLBLFile(lblname)
  writeLog(_log, 'Finished reading label file...')
  #
  # Parse important data from label file
  #
  InstrPresum = lblDic['INSTR_MODE_ID']['Presum']	# Onboard presums
  fft_length = fft_length / InstrPresum
  instrMode = lblDic['INSTR_MODE_ID']['Mode']
  BitsPerSample = lblDic['INSTR_MODE_ID']['BitsPerSample']
  nrec = lblDic['FILE_RECORDS']				# Number of records
  #
  # Read in Auxilliary file data 
  # I should probably read this in chunks rather than the whole file at once
  #
  writeLog(_log, 'Reading Auxilliary file...')
  AuxDF = parseAuxFile(auxname)  
  writeLog(_log, 'Finished reading Auxilliary file...')
  ###############################################################################################
  # 
  # End input data
  # No other external variables required
  #
  ############################### BEGIN PROCESSING ##############################################
  if InstrPresum == 1:
    presum_proc = 4.0					# Integrate 4 pulses at a time to limit data volume
    fft_length = 64					# Effective length of 256 pulses
    nsums = nrec / 4.0					# Number of blocks of 4 records to read
  else:							# If instrument presum not 1
    presum_proc = 1.0					# Keep presum values the same as on the instrument
    fft_length = 256.0 / InstrPresum			# Keep a total of 256 pulses in any integration (uniform Doppler bin size)
    nsums = np.copy(nrec)
  data = np.zeros([3600,1], float)			# Partial data array (columns = pulse record; rows = time delay)
  presum_total = InstrPresum * presum_proc		# Final compression factor of dechirped data
  new_x = np.zeros(3600, float)				# For lining up traces
  #############################################################################
  #
  if verb:
    writeLog(_log, 'Number of Records:\t{}'.format(nrec))
    writeLog(_log, 'InstrPresum:\t{}'.format(InstrPresum))
    writeLog(_log, 'Instrument Mode:\t{}'.format(instrMode))
  #############################################################################
  #
  # Construct the window function
  #
  window, win_str = makeWindow(win_type)
  writeLog(_log, win_str)
  #############################################################################
  #
  # Determine the record length based on the bit resolution of EDR data
  #
  if BitsPerSample == 4:
    recLen = 1986
  elif BitsPerSample == 6:
    recLen = 2886
  elif BitsPerSample == 8:
    recLen = 3786
  tmp = np.zeros(1800, complex)
  _file1 = open(edrname, 'rb')					# Open EDR file
  Fc = (80/3 - 20) 						# MHz
  dt = 3/80							# Microseconds
  t = np.arange(0,4096*dt, dt)					# Time vector
  C = np.exp(2*np.pi*1j*Fc*t)					# Complex exponential
  decomp = np.zeros([2048, nsums], complex)			# Matrix for range compressed data
  timeArr = np.zeros(nsums, float)				# Time at each pulse
  latArr = np.zeros(nsums, float)				# Latitude array 
  lonArr = np.zeros(nsums, float) 				# Longitude array
  radiusArr = np.zeros(nsums, float)				# Radius value at each pulse
  tanVelArr = np.zeros(nsums, float)				# tangential velocity
  radVelArr = np.zeros(nsums, float)				# radial velocity
  position = np.zeros([nsums, 3], float)			# position matrix
  velocity = np.zeros([nsums, 3], float)			# velocity matrix
  ##################################### BEGIN RANGE PROCESSING ############################################
  # 
  # If data is presum-1
  #
  if InstrPresum == 1:					# Do a running 4X presum
    if verb:
      writeLog(_log, 'Begin decompression at:{}'.format(datatime.now()))
      for _i in range(0,nsums):
        #
        # Determine and load calibrated chirp
        #
        calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i])
        timeArr = np.zeros(nsums, float)				# Time at each pulse
        latArr = np.zeros(nsums, float)					# Latitude array 
        lonArr = np.zeros(nsums, float) 				# Longitude array
        radiusArr = np.zeros(nsums, float)				# Radius value at each pulse
        tanVelArr = np.zeros(nsums, float)				# tangential velocity
        radVelArr = np.zeros(nsums, float)				# radial velocity
        position = np.zeros([nsums, 3], float)				# position matrix
        velocity = np.zeros([nsums, 3], float)				# velocity matrix
        for _k in range(0,4):
          #
          # Read in a record
          #
          sci, ancil = readEDRrecord(_file1, _k+4*_i, recLen, BitsPerSample) 
          ancil = parseAncilliary(ancil)
          #
          # Fill in necessary arrays for averaging
          #
#          t[_k] = ancil['TIME_N']
          t[_k] = AuxDF['ELAPSED_TIME'][4*_i+_k]
          lat[_k] = AuxDF['SUB_SC_PLANETOCENTRIC_LATITUDE'][4*_i+_k]
          lon[_k] = AuxDF['SUB_SC_EAST_LONGITUDE'][4*_i+_k]
          r[_k] = ancil['RADIUS_N']
          vt[_k] = ancil['TANGENTIAL_VELOCITY_N']
          vr[_k] = ancil['RADIAL_VELOCITY_N']
          p[_k,:] = [ AuxDF['X_MARS_SC_POSITION_VECTOR'][4*_i+_k],
                      AuxDF['Y_MARS_SC_POSITION_VECTOR'][4*_i+_k],
                      AuxDF['Z_MARS_SC_POSITION_VECTOR'][4*_i+_k]
                    ]
          v[_k,:] = [ AuxDF['X_MARS_SC_VELOCITY_VECTOR'][4*_i+_k],
                      AuxDF['Y_MARS_SC_VELOCITY_VECTOR'][4*_i+_k],
                      AuxDF['Z_MARS_SC_VELOCITY_VECTOR'][4*_i+_k]
                    ]
          #
          # Prepare data for range compression
          #
          echoes = np.zeros(4096, float)
          echoes[:len(sci)] = sci				# Pad with zeros
          echoes = echoes * C					# Multiply by complex exponential
          ecSpec = np.fft.fft(echoes)				# Compute the FFT
          ecSpec = np.fft.fftshift(ecSpec)			# Shift the data in ascending frequency order
          ecSpec = ecSpec[1024:3072]				# Select only the central 2048 samples
          if fil_type == 'Match':
              decomp[:,_i] += np.fft.ifft(np.fft.ifftshit(window*calChirp * ecSpec))
          elif fil_type == 'Inverse':
              decomp[:,-i] += np.fft.ifft(np.fft.ifftshift(window*(ecSpec / calChirp)))
        #
        # Fill in necessary data from either the auxilliary or ancilliary data
        # for each pulse
        #
        timeArr[_i] = np.average(t)				# Time at each pulse
        latArr[_i] = np.average(lat)				# Latitude array 
        lonArr[_i] = np.average(lon) 				# Longitude array
        radiusArr[_i] = np.average(r)				# Radius value at each pulse
        tanVelArr[_i] = np.average(vt)				# tangential velocity
        radVelArr[_i] = np.average(vr)				# radial velocity
        position[_i,:] = [ np.average(p[:,0]),			# position arrays
                     np.average(p[:,1]),
                     np.average(p[:,2])
                   ]			
        velocity = [ np.average(v[:,0]),			# velocity arrays
                     np.average(v[:,1]),
                     np.average(v[:,2])
                   ]
  ##################################### END PRESUM-1 PROCESSING ###########################################
  if InstrPresum != 1:				
    if verb:
      writeLog(_log, 'Begin decompression at:\t{}'.format(datetime.now()))
    for _i in range(0, nsums):
      sci, ancil = readEDRrecord(_file1, _i, recLen, BitsPerSample)  
      ancil = parseAncilliary(ancil)
      #
      # Fill in necessary data from either the auxilliary or ancilliary data
      # for each pulse
      #
#      timeArr[_i] = ancil['TIME_N']
      timeArr[_i] = AuxDF['ELAPSED_TIME'][_i]
      latArr[_i] = AuxDF['SUB_SC_PLANETOCENTRIC_LATITUDE'][_i]
      lonArr[_i] = AuxDF['SUB_SC_EAST_LONGITUDE'][_i]
      radiusArr[_i] = ancil['RADIUS_N']
      tanVelArr[_i] = ancil['TANGENTIAL_VELOCITY_N']
      radVelArr[_i] = ancil['RADIAL_VELOCITY_N']
      position[_i,:] = [ AuxDF['X_MARS_SC_POSITION_VECTOR'][_i],
                         AuxDF['Y_MARS_SC_POSITION_VECTOR'][_i],
                         AuxDF['Z_MARS_SC_POSITION_VECTOR'][_i]
                        ]			
      velocity[_i,:] = [ AuxDF['X_MARS_SC_VELOCITY_VECTOR'][_i],
                         AuxDF['Y_MARS_SC_VELOCITY_VECTOR'][_i],
                         AuxDF['Z_MARS_SC_VELOCITY_VECTOR'][_i]
                       ]
      #
      # Determine and load calibrated chirp
      #
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i])
      sys.exit()
      #
      # Prepare data for range compression
      # 
      echoes = np.zeros(4096, float)
      echoes[:len(sci)] = sci					# Pad with zeros
      echoes = echoes * C					# Multiply by complex exponential
      ecSpec = np.fft.fft(echoes)				# Compute the FFT
      ecSpec = np.fft.fftshift(ecSpec)				# Shift the data in ascending frequency order
      ecSpec = ecSpec[2048:] + 1j*ecSpec[:2048]			# Thought this was worth a shot
#      ecSpec = ecSpec[1024:3072]				# Select only the central 2048 samples
      if fil_type == "Match":
        decomp[:, _i] = np.fft.ifft(np.fft.ifftshift(window*(calChirp * ecSpec)))
      elif fil_type == "Inverse":
        decomp[:, _i] = np.fft.ifft(np.fft.ifftshift(window*(ecSpec / calChirp)))
  ################################ END OF PULSE COMPRESSION AND DECHIRP OPERATIONS ########################
  writeLog(_log, 'Decompression ended at:\t{}'.format(datetime.now()))
  ################################ Begin Azimuth Processing ##############################################
  #
  # Determine number of output rows
  #
  np.save('decomp.npy', decomp)
#  decomp = np.load('decomp.npy')
#  plt.imshow(np.real(decomp))
#  plt.show()
  dtor = np.pi/180.
  central_angle = 2.0*np.arcsin(np.sqrt(np.sin((np.abs(latArr[0] - latArr[-1]))*dtor/2.0)**2 + \
                  np.cos(latArr[0]*dtor)*np.cos(latArr[-1]*dtor)* \
                  np.sin((np.abs(lonArr[0]-lonArr[-1]))*dtor/2.0)**2))
  num_output_rows = int((central_angle/dtor)/0.0078125)
  complex_out = np.empty([3600, num_output_rows], complex)
#  map_out = np.empty([3600, num_output_rows, 3], float)				# RGB File
  pow_out = np.empty([3600, num_output_rows], float)					# BW Power File
  temp1 = np.zeros([2048, int(fft_length)], complex)					# Chunk of voltage array
  if verb:
    writeLog(_log, 'Begin summation to output arrays at:\t{}'.format(datetime.now()))
  for _k in range(0, num_output_rows):
    time_row = timeArr[_k]					# This is not correct
    #
    # Define Doppler Window of length FFT_LENGTH 
    #
    dopp_bin_center = int((0.5 + time_row * (700.28 / presum_total)))
    min_dopp = dopp_bin_center - fft_length/2.0
    max_dopp = dopp_bin_center + fft_length/2.0
    if max_dopp >= nsums:
      continue							# Off right edge
    if min_dopp >= 0:
      temp1[:,:] = decomp[:,:fft_length]

  _log.close()
#  plt.imshow(np.real(decomp), cmap='gray')
#  plt.show()
  return

if __name__ == '__main__':
  main()
