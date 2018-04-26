import numpy as np

def rangeCompression(sci, calChirp, window, fil_type='Match'):
  """
    This function performs the chirp compression on a single SHARAD record
    Inputs:
      sci - range (fast-time) record of decompressed science data
      calChirp - Complex conjugate of the cal. chirp spectrum
      fil_type  - Match or Inverse
  """
  Fc = (80/3 - 20)                                      # MHz
  dt = 3/80                                             # Microseconds
  t = np.arange(0,4096*dt, dt)                          # Time vector
  C = np.exp(2*np.pi*1j*Fc*t)                           # Complex exponential
  #
  # Check length of the science data
  #
  if len(sci) != 4096:
    echoes = np.zeros(4096, complex)
    echoes[:len(sci)] = sci
  echoes = echoes #* C
  #
  # Compute the FFT
  #
  ecSpec = np.fft.fft(echoes)
  #
  # Shift the data into ascending frequency order
  #
  ecSpec = np.fft.fftshift(ecSpec)
  #
  # Take central 2048 samples
  #
  ecSpec = ecSpec[1024:3072]
  #
  # Perform Chirp compression
  #
  if fil_type == 'Match':
    temp = window * calChirp * ecSpec
  elif fil_type == "Inverse":
    temp = window * (ecSpec / calChirp)
  temp2 = np.fft.ifftshift(temp)
  decomp = np.fft.ifft(temp2)
  return decomp
