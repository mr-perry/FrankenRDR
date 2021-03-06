import numpy as np
import matplotlib.pyplot as plt

def rangeCompression(sci, calChirp, window, chirp='ref', fil_type='Match', diag=False):
  """
    This function performs the chirp compression on a single SHARAD record
    Inputs:
      sci - range (fast-time) record of decompressed science data
      calChirp - Cal. Chirp spectrum
      fil_type  - Match or Inverse
  """
  #
  # Get Chirp Length to allow use of the ideal chirp
  # 2048 -- reference chirp
  # 3600 -- ideal chirp
  #
  length = len(calChirp)
  if chirp == 'ideal' or chirp == 'UPB':
    dt = 135.00e-6 / 3600.
    sharad_ipp = 1.0 / 700.28
    Fc = 1 / dt
    plen = 85.05e-6
    nsamp = plen / dt
    t = np.arange(0,length) * dt
    ecSpec = np.fft.fft(sci) / len(sci)
    ecFreq = np.fft.fftshift(np.fft.fftfreq(length, dt)) 
    ecSpec_shift = np.fft.fftshift(ecSpec)
    ecFreq_shift = np.fft.fftshift(ecFreq) 
    dechirp = np.conj(calChirp)
    decomp =  np.fft.ifft(window*(dechirp*(ecSpec))) * len(sci)
    if diag:
      #
      # Plot original time signal
      #
      plt.subplot(3,1,1)
      plt.plot(t*1e6, np.abs(np.real(sci)))
      plt.title("Raw Signal")
      plt.xlabel('Time (s)')
      plt.ylabel('Amplitude')
      #
      # Spectra
      # 
      plt.subplot(3,1,2)
      plt.plot(ecFreq/1e6, np.abs(np.real(ecSpec_shift)))
  #    plt.plot(ecFreq, window)
      plt.title("Windowed Magnitude Spectrum")
      plt.subplot(3,1,3)
      plt.plot(ecFreq/1e6, np.unwrap(np.angle(ecSpec_shift)))
      plt.title('Angle Spectrum')
      plt.tight_layout()
      plt.show()
      plt.clf()
      #
      # Chirp Spectra
      #
      chirpFreq = np.linspace(-(13+1/3), (13+1/3),num=3600)
      plt.subplot(2,1,1)
      plt.plot(chirpFreq, np.abs(np.real(calChirp)))
      plt.title('Chirp Spectrum')
      plt.subplot(2,1,2)
      plt.plot(chirpFreq, np.unwrap(np.angle(calChirp)))
      plt.title('Chirp Phase Spectrum')
      plt.show()
      plt.clf()
      #
      # Range Compression
      #
      plt.subplot(3,1,1)
      plt.plot(ecFreq, np.abs(np.real((dechirp))))
      plt.title('Dechirp Power Spectrum')
      plt.subplot(3,1,2)
      plt.plot(ecFreq, np.unwrap(np.angle(dechirp))) 
      plt.title('Dechirp Phase Spectrum')
      plt.subplot(3,1,3)
      plt.plot(t, np.abs(np.real(decomp)))
      plt.title('Range Compressed Amplitude')
      plt.tight_layout()
      plt.savefig('diag_plot.png', dpi=1000)
      plt.show()
      sys.exit()
    return decomp
  elif chirp == 'ref' or chirp == 'vib':
    dt = (3./80.)                                           # 0.0375 Microseconds
    t = np.arange(0*dt, 4096*dt, dt)
    #
    # Check length of the science data
    #
    if len(sci) != 4096:
      echoes = np.zeros(4096, complex)
      echoes[:len(sci)] = sci
    if chirp == 'ref':
      #
      # Check the length of the window
      #
      if len(window) != 4096:
        temp = np.copy(window)
        window = np.zeros(4096, float)
        window[:len(temp)] = temp
    elif chirp == 'vib': # no window
      window = np.ones(4096, float)
    #
    # Compute the FFT and place on same scale as the calibrated chirp
    #
    ecSpec = np.fft.fft(echoes)  / len(echoes)
    ecFreq = np.fft.fftfreq(len(echoes), d=dt)
    if chirp == 'ref':
      #
      # Pad Chirp
      # 
      chirpWin = np.zeros(4096, complex)
      chirpWin[0:2048] = np.conj(calChirp)
    elif chirp == 'vib':
      chirpWin = np.conj(calChirp)
    #
    # Perform Chirp compression
    #
    if fil_type == 'Match':
      dechirp = window * ecSpec * chirpWin
    elif fil_type == "Inverse":
      sys.exit()
    #
    # Inverse Fourier transform and fix scaling
    #
    decomp = np.fft.ifft(dechirp) * len(dechirp)
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
      plt.plot(np.abs(ecSpec))
      plt.title("Modulus Amplitude Spectrum (|S|)")
      plt.subplot(3,1,3)
      plt.plot(np.unwrap(np.angle(ecSpec)))
      plt.title('Unwrapped Phase Spectrum')
      plt.tight_layout()
      plt.show()
      plt.clf()
      #
      # Window
      #
      plt.plot(window)
      plt.title('Windowing Function')
      plt.show()
      plt.clf()
      #
      # Chirp Spectra
      #
      plt.subplot(2,1,1)
      plt.plot(np.abs(chirpWin))
      plt.title('Modulus Amplitude Spectrum of Reference Chirp (|C|)')
      plt.subplot(2,1,2)
      plt.plot(np.unwrap(np.angle(chirpWin)))
      plt.title('Unwrapped Phase Spectrum of Reference Chirp')
      plt.tight_layout()
      plt.show()
      plt.clf()
      #
      # Range Compression
      #
      plt.subplot(3,1,1)
      plt.plot(np.abs(dechirp))
      plt.title('Modulus Amplitude Spectrum of Dechirped Data (| S * C |')
      plt.subplot(3,1,2)
      plt.plot(np.unwrap(np.angle(dechirp))) 
      plt.title('Unwrapped Phase Spectrum of Dechirped Data')
      plt.subplot(3,1,3)
      plt.plot(np.real(decomp))
      plt.title('Time Series Amplitude of Dechirped Data IFFT(| S * C |)')
      plt.tight_layout()
      plt.show()
      #
      # Other plots
      #
      plt.subplot(3,1,1)
      plt.plot(np.real(decomp))
      plt.title('Real Component of Dechirped Time Series')
      plt.subplot(3,1,2)
      plt.plot(np.imag(decomp))
      plt.title('Imaginary Component of Dechirped Time Series')
      plt.subplot(3,1,3)
      plt.plot(np.unwrap(np.angle(decomp)))
      plt.title('Unwrapped Phase of Dechirped Time Series')
      plt.tight_layout()
      plt.show()
    return decomp
