import numpy as np
import matplotlib.pyplot as plt

def rangeCompression(sci, calChirp, window, fil_type='Match', diag=False):
  """
    This function performs the chirp compression on a single SHARAD record
    Inputs:
      sci - range (fast-time) record of decompressed science data
      calChirp - Complex conjugate of the cal. chirp spectrum
      fil_type  - Match or Inverse
  """
  Fc = (80./3. - 20.)                                      # MHz
  dt = (3./80.)                                           # Microseconds
  t = np.arange(0,4096*dt, dt)                          # Time vector
  arg = np.arange(0,4096) * 2.0 * np.pi * Fc*1.0e6 * dt*1e-6
  pha_shift = np.cos(arg) + 1j*np.sin(arg)		 # Complex exponential using Euler's formula
  #
  # Check length of the science data
  #
  if len(sci) != 4096:
    echoes = np.zeros(4096, complex)
    echoes[:len(sci)] = sci
  echoes = echoes * pha_shift
  #
  # Compute the FFT
  #
  ecSpec = np.fft.fft(echoes)
  ecFreq = np.fft.fftfreq(len(echoes), d=dt*10**-6)
  #
  # Take central 2048 samples
  #
  st = 1024
  en = 3072
  ecSpec = ecSpec[st:en]
  ecFreq = ecFreq[st:en]
  #
  # Multiply chirp by length of 2048 to put in same FFT units as Python FFT
  #
  calChirp = calChirp * len(ecSpec)
  #
  # Perform Chirp compression
  #
  if fil_type == 'Match':
    temp = window * np.conj(calChirp) * ecSpec
  elif fil_type == "Inverse":
    temp = window * (ecSpec / np.conj(calChirp))
  #
  # Inverse Fourier transform
  #
  decomp = np.fft.ifft(temp) 
  if diag:
    plt.subplot(4,1,1)
    plt.plot(np.power(np.abs(ecSpec), 2))
    plt.subplot(4,1,2)
    plt.plot(np.power(np.abs(calChirp), 2))
    plt.subplot(4,1,3)
    plt.plot(np.power(np.abs(temp),2 ))
    plt.subplot(4,1,4)
    plt.plot(np.power(np.abs(decomp),2 ))
    plt.show()
    sys.exit()
  return decomp
