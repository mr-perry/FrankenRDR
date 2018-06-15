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
  pha_shift = np.exp(2*np.pi*1j*Fc*t)
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
  ecSpec = np.fft.fftshift(ecSpec)
  ecFreq = np.fft.fftshift(ecFreq)
  ecSpec_win = window * ecSpec
  #
  # Take central 2048 samples
  #
  st = 1024
  en = 3072
  if diag:
    print(ecFreq[st], ecFreq[en])
  ecSpec_cut = ecSpec_win[st:en]
  ecFreq_cut = ecFreq[st:en]
  #
  # Perform Chirp compression
  #
  if fil_type == 'Match':
#    temp = window * np.conj(calChirp) * ecSpec_cut
    dechirp = (ecSpec_cut) * calChirp 
#    dechirp = (ecSpec_cut) * np.conj(calChirp) 
  elif fil_type == "Inverse":
    sys.exit()
  #
  # Pad spectrum with zeros
  #
  dechirp_pad = np.zeros(4096, complex)
  dechirp_pad[1024:3072] = dechirp
  #dechirp_pad = dechirp
  #
  # Shift the data back to standard order
  #
  shift_dechirp_pad = np.fft.ifftshift(dechirp_pad)
  #
  # Inverse Fourier transform and fix scaling
  #
  decomp = np.fft.ifft(shift_dechirp_pad) * len(shift_dechirp_pad) 
  if diag:
    #
    # Plot original time signal
    #
    plt.subplot(5,1,1)
    plt.plot(t, np.abs(np.real(echoes)))
    plt.title("Raw Signal")
    plt.xlabel('Time (s)')
    plt.ylabel('"Power"')
    #
    # Spectra
    # 
    plt.subplot(5,2,3)
    plt.plot(ecFreq, np.abs(np.real(ecSpec)))
#    plt.xlim(-(6+2/3),(6+2/3))
    plt.plot(ecFreq, window)
    plt.title("Windowed Magnitude Spectrum")
    plt.subplot(5,2,4)
    plt.plot(ecFreq, np.angle(ecSpec))
    plt.xlim(-(6+2/3),(6+2/3))
    plt.title('Angle Spectrum')
    #
    # Chirp Spectra
    #
    chirpFreq = np.arange(0*0.00651042*2-(13+1/3), 2048*0.00651042*2-(13+1/3), 0.00651042*2)
    plt.subplot(5,2,5)
    plt.plot(chirpFreq, np.abs(np.real(calChirp)))
    plt.title('Chirp Spectrum')
    plt.subplot(5,2,6)
    plt.plot(chirpFreq, np.angle(calChirp))
    plt.title('Chirp Phase Spectrum')
    #
    # Range Compression
    #
    plt.subplot(5,2,7)
    plt.plot(ecFreq_cut, np.abs(np.real((dechirp))))
    plt.title('Dechirp Power Spectrum')
    plt.subplot(5,2,8)
    plt.plot(ecFreq_cut, np.angle(dechirp)) 
    plt.title('Dechirp Phase Spectrum')
    plt.subplot(5,1,5)
    plt.plot(np.abs(np.real(decomp)))
    plt.title('Range Compressed Amplitude')
    plt.tight_layout()
    plt.savefig('diag_plot.png', dpi=1000)
    plt.show()
    sys.exit()
  return decomp
