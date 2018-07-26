#
# Import necessary libraries
#
import numpy as np
import os, sys
import matplotlib.pyplot as plt


def detChirpFiles(TxTemp, RxTemp, chirp='ref'):
  """
  This function determines the appropriate calibrated chirp file to use
  for range compression and returns the decoded calibrated chirp.
  This is solely based off the temperatures of the TxTemp and RxTemp.
  """
  if chirp == 'ref' or chirp == 'vib':
    calibRoot = '../calib/'
    calibName = 'reference_chirp'
    ext = '.dat'
    TxCalNames = ['m20tx', 'm15tx', 'm10tx', 'm05tx',
                  'p00tx', 'p20tx', 'p40tx', 'p60tx']
    RxCalNames = ['m20rx', 'p00rx', 'p20rx', 'p40rx',
                  'p60rx']
    #
    # Define vectors for Tx and Rx temps
    #
    Tx = [-20, -15, -10, -5, 0, 20, 40, 60]
    Rx = [-20, 0, 20, 40, 60]
    calibChirpFiles = []
    TxDiff = []
    RxDiff = []
    #
    # Find distance
    #
    TxDiff[:] = [abs(x - TxTemp) for x in Tx]
    RxDiff[:] = [abs(x - RxTemp) for x in Rx]
    #
    # Find the indices of the closest Tx and Rx value
    #
    calibTx = TxCalNames[TxDiff.index(min(TxDiff))]
    calibRx = RxCalNames[RxDiff.index(min(RxDiff))]
    #
    # Construct File name
    #
    calChirpFile = calibRoot + calibName + '_' + \
                   TxCalNames[TxDiff.index(min(TxDiff))] + '_' + \
                   RxCalNames[RxDiff.index(min(RxDiff))] + ext
    if os.path.isfile(calChirpFile):
      calChirp = np.fromfile(calChirpFile, dtype='<f')
      if chirp == 'ref':
        real = calChirp[:2048]
        imag = calChirp[2048:]
      elif chirp == 'vib':
        #
        # Add a zero to the end to mimic missing Nyquist
        #
        real = np.zeros(4096, float)
        real[0:2048] = calChirp[:2048]
        #
        # Drop the last sample to properly stitch Nyquist and reverse
        #
        real[2049:] = np.flipud(calChirp[1:2048])
        #
        # Add a zero to the front
        #
        imag = np.zeros(4096, float)
        imag[0:2048] = calChirp[2048:]
        #
        # Drop the last sample and reverse and change sign
        #
        imag[2049:] = -1 * np.flipud(calChirp[2049:])
      calChirp = real + 1j*imag
      return calChirp
    else:
      print('Calibrated chirp file not found...exiting.')
      sys.exit()
  else:
    #
    # Use ideal chirp from UPB
    #
    delay_res = 135.00e-6 / 3600.000
    sharad_ipp = 1.0 / 700.28
    fhi = 25.00e6
    flo = 15.00e6
    plen = 85.05e-6
    nsamp = plen / delay_res
    fslope = (flo - fhi) / plen
    ctime = np.arange(0,nsamp) * delay_res
    arg = 2.0*np.pi*ctime*(fhi+fslope*ctime/2.0)
    ideal_chirp = np.zeros(3600, complex)
    ideal_chirp[:int(nsamp)] = np.sin(arg)
    ideal_chirp_FFT = np.fft.fft(ideal_chirp)
    if chirp == 'UPB':
      #
      # Load cal_filter.dat
      #
      cal_filter = np.fromfile('../calib/cal_filter.dat', '<f')
      cal_filter = cal_filter[:1800] + 1j*cal_filter[1800:]
      cal_filter = np.roll(cal_filter, 900)
      calChirp = np.zeros(3600, complex)
      calChirp[1800:] = cal_filter
      calChirp = calChirp*ideal_chirp_FFT
      return calChirp
    else:
      calChirp = ideal_chirp_FFT
      return calChirp
