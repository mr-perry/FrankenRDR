import numpy as np
import matplotlib.pyplot as plt

def rangeCompression(sci, calChirp, window, fil_type='Match', diag=False):
  """
    This function performs the chirp compression on a single SHARAD record
    Inputs:
      sci - range (fast-time) record of decompressed science data
      calChirp - Cal. Chirp spectrum
      fil_type  - Match or Inverse
  """
  Fc = (80./3. - 20.)                                    # MHz
  dt = (3./80.)                                           # Microseconds
  t = np.arange(0*dt, 4096*dt, dt)
  pha_shift = np.exp(2*np.pi*complex(0,1)*Fc*t)
  #
  # Check length of the science data
  #
  if len(sci) != 4096:
    echoes = np.zeros(4096, complex)
    echoes[:len(sci)] = sci
  echoes_shift = echoes * pha_shift
  #
  # Compute the FFT and place on same scale as the calibrated chirp
  #
  ecSpec = np.fft.fft(echoes_shift) / len(echoes_shift)
  ecFreq = np.fft.fftfreq(len(echoes_shift), d=dt)
  #
  # Put the Spectra into natural order
  #
  ecSpec = np.fft.fftshift(ecSpec)
  ecFreq = np.fft.fftshift(ecFreq)
  #
  # Take central 2048 samples
  #
  st = 1024
  en = 3072
  ecSpec_cut = ecSpec[st:en]
  ecFreq_cut = ecFreq[st:en]
  #
  # Perform Chirp compression
  #
  if fil_type == 'Match':
    dechirp = (ecSpec_cut * window) * (calChirp)
    #dechirp = (ecSpec_cut * window) * np.conj(calChirp)
  elif fil_type == "Inverse":
    sys.exit()
  #
  # Shift the data back to standard order
  #
  shift_dechirp = np.fft.ifftshift(dechirp)
  #
  # Inverse Fourier transform and fix scaling
  #
  decomp = np.fft.ifft(shift_dechirp) * len(shift_dechirp) 
  if diag:
    #
    # Plot original time signal
    #
    plt.subplot(3,1,1)
    plt.plot(t, np.abs(np.real(echoes)))
    plt.title("Raw Signal")
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    #
    # Spectra
    # 
    plt.subplot(3,1,2)
    plt.plot(ecFreq, np.abs(np.real(ecSpec)))
#    plt.xlim(-(6+2/3),(6+2/3))
#    plt.plot(ecFreq, window)
    plt.title("Windowed Magnitude Spectrum")
    plt.subplot(3,1,3)
    plt.plot(ecFreq, np.unwrap(np.angle(ecSpec)))
    plt.xlim(-(6+2/3),(6+2/3))
    plt.title('Angle Spectrum')
    plt.show()
    plt.clf()
    #
    # Chirp Spectra
    #
    chirpFreq = np.linspace(-(13+1/3), (13+1/3),num=4096 )
    plt.subplot(2,1,1)
    plt.plot(chirpFreq[1024:3072], np.abs(np.real(calChirp)))
    plt.title('Chirp Spectrum')
    plt.subplot(2,1,2)
    plt.plot(chirpFreq[1024:3072], np.unwrap(np.angle(calChirp)))
    plt.title('Chirp Phase Spectrum')
    plt.show()
    plt.clf()
    #
    # Range Compression
    #
    plt.subplot(3,1,1)
    plt.plot(ecFreq_cut, np.abs(np.real((dechirp))))
    plt.title('Dechirp Power Spectrum')
    plt.subplot(3,1,2)
    plt.plot(ecFreq_cut, np.unwrap(np.angle(dechirp))) 
    plt.title('Dechirp Phase Spectrum')
    plt.subplot(3,1,3)
    plt.plot(np.abs(np.real(decomp)))
    plt.title('Range Compressed Amplitude')
    plt.tight_layout()
    plt.savefig('diag_plot.png', dpi=1000)
    plt.show()
    sys.exit()
  return decomp
