import numpy as np
import matplotlib.pyplot as plt


def main():
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
  ideal_chirp[:len(arg)] = np.sin(arg)
  ideal_chirp_FFT = np.fft.fft(ideal_chirp)
  ideal_chirp_freqs = np.fft.fftfreq(3600, delay_res)
  #
  # Load nal_filter.dat
  #
  cal_filter = np.fromfile('../calib/cal_filter.dat', '<f')
  cal_filter = cal_filter[:1800] + 1j*cal_filter[1800:]
  cal_filter = np.roll(cal_filter, 900)
  dechirp = np.zeros(3600, complex)
  dechirp[1800:] = cal_filter
  dechirp = (dechirp*ideal_chirp_FFT)
  dechirp_IFFT = np.fft.ifft(dechirp)
  #
  # Load a reference chirp example
  #
  raw_ref_chirp = np.fromfile('../calib/reference_chirp_p00tx_p00rx.dat', '<f')
  #ref_chirp = np.conj(ref_chirp[:2048] + 1j*ref_chirp[2048:])
  #raw_ref_chirp = (raw_ref_chirp[:2048] + 1j*raw_ref_chirp[2048:])
  raw_real = raw_ref_chirp[2047::-1]
  raw_imag = raw_ref_chirp[4096:2047:-1]
  raw_ref_chirp = raw_real + 1j*raw_imag
  ref_chirp = np.zeros(4096, complex)
  ref_chirp[2048:] = raw_ref_chirp
  plt.subplot(2,2,1)
  plt.plot(np.real(ref_chirp[2048:]))
  plt.title('Reference Chirp -- Real')
  plt.subplot(2,2,2)
  plt.plot(np.imag(ref_chirp[2048:]))
  plt.title('Reference Chirp -- Imaginary')
  plt.subplot(2,2,3)
  plt.plot(np.real(dechirp[1552:]))
  plt.title('UPB Dechirp -- Real')
  plt.subplot(2,2,4)
  plt.plot(np.imag(dechirp[1552:]))
  plt.title('UPB Dechirp -- Imaginary')
  plt.tight_layout()
  plt.show()
  #
  # Let's plot
  #
  plt.subplot(3,2,1)
  plt.plot(np.real(np.fft.ifft(np.fft.ifftshift(dechirp))))
  plt.title('Dechirp - Time Series Real')
  plt.subplot(3,2,3)
  plt.plot(np.power(np.abs(dechirp),2))
  plt.title('Dechirp - Power Spectrum')
  plt.subplot(3,2,5)
  plt.plot(np.unwrap(np.angle(dechirp)))
  plt.title('Dechirp - Phase Spectrum')
  #
  # Let's plot
  #
  plt.subplot(3,2,2)
  plt.plot((np.real(np.fft.ifft(np.fft.ifftshift(ref_chirp)))))
  plt.title('Reference Chirp - Time Series Real')
  plt.subplot(3,2,4)
  plt.plot(np.power(np.abs(ref_chirp),2))
  plt.title('Reference Chirp - Power Spectrum')
  plt.subplot(3,2,6)
  plt.plot(np.unwrap(np.angle(ref_chirp)))
  plt.title('Reference Chirp - Phase Spectrum')
  plt.tight_layout()
  plt.show()
  #
  # Plot real and imaginary parts of the PDS chirp and cal_filter
  #
  '''
  plt.subplot(2,2,1)
  plt.plot(np.abs(np.real(raw_ref_chirp)))
  plt.title('Raw Reference Chirp - Amplitude')
  plt.subplot(2,2,2)
  plt.plot(np.unwrap(np.angle(raw_ref_chirp)))
  plt.title('Raw Reference Chirp - Phase (Unwrapped)')
  plt.subplot(2,2,3)
  plt.plot(np.abs(np.real(dechirp[1552:])))
  plt.title('cal_filter * fft(ideal_chirp) -- Amplitude')
  plt.subplot(2,2,4)
  plt.plot(np.unwrap(np.angle(dechirp[1552:])))
  plt.title('cal_filter * fft(ideal_chirp) -- Phase (Unwrapped)')
  plt.tight_layout()
  plt.show()
  '''
  return

if __name__ == '__main__':
  main()
