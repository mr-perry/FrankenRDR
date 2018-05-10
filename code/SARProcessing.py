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
  Fc = (80./3. - 20.)*1.e6                                    # MHz
  dt = (3./80.)*1.e-6                                           # Microseconds
  t = np.arange(0,4096*dt, dt)                          # Time vector
  arg = np.arange(0,4096) * 2.0 * np.pi * Fc * dt
  pha_shift = np.cos(arg) + 1j*np.sin(arg)		 # Complex exponential using Euler's formula
  #
  # Check length of the science data
  #
  if len(sci) != 4096:
    echoes = np.zeros(4096, complex)
    echoes[:len(sci)] = sci
  echoes = echoes * pha_shift
  #
  # Compute the FFT and place on same scale as the calibrated chirp
  #
  ecSpec = np.fft.fft(echoes) / len(echoes)
  ecFreq = np.fft.fftfreq(len(echoes), d=dt)
  if diag:
    print(ecFreq[1024], ecFreq[2047], ecFreq[2048], ecFreq[3072])
  #
  # Take central 2048 samples
  #
  st = 1024
  en = 3072
  ecSpec = ecSpec[st:en]
  ecFreq = ecFreq[st:en]
  #
  # Shift the window to maximum of ecSpec
  # 
  #mx = np.where(ecSpec == np.amax(ecSpec))
  #md_wn = len(window)/2 - 1
  #sft = int(mx[0][0] - md_wn)
  #window = np.roll(window, sft) 
  #
  # Perform Chirp compression
  #
  if fil_type == 'Match':
    temp = window * np.conj(calChirp) * ecSpec
  elif fil_type == "Inverse":
    temp = window * (ecSpec / np.conj(calChirp))
  #
  # Inverse Fourier transform and fix scaling
  #
  decomp = np.fft.ifft(temp) * len(temp) 
  if diag:
    plt.subplot(5,1,1)
    plt.plot(np.power(np.abs(ecSpec), 2))
    plt.subplot(5,1,2)
    plt.plot(window)
    plt.subplot(5,1,3)
    plt.plot(np.power(np.abs(calChirp), 2))
    plt.subplot(5,1,4)
    plt.plot(np.power(np.abs(temp),2 ))
    plt.subplot(5,1,5)
    plt.plot(np.power(np.abs(decomp),2 ))
    plt.show()
  return decomp
