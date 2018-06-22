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
#  plt.subplot(2,1,1)
#  plt.plot((np.real(ideal_chirp)))
#  plt.title('Ideal Chirp Time Series - Real')
#  plt.subplot(2,1,2)
#  plt.plot((np.imag(ideal_chirp)))
#  plt.title('Ideal Chirp Time Series - Imag')
#  plt.xlabel('Sample')
#  plt.tight_layout()
#  plt.savefig('../runs/FritzFigsRequest/IdealChirp_TimeSeries.png', dpi=1000)
  ideal_chirp_FFT = np.fft.fft(ideal_chirp)
  ideal_chirp_freqs = np.fft.fftfreq(3600, delay_res)
#  plt.clf()
#  plt.subplot(4,1,1)
#  plt.plot(np.real(ideal_chirp_FFT))
#  plt.title('Ideal Chirp Spectra - Real (Standard Order)')
#  plt.subplot(4,1,2)
#  plt.plot(np.imag(ideal_chirp_FFT))
#  plt.title('Ideal Chirp Spectra - Imaginary (Standard Order)')
#  plt.subplot(4,1,3)
#  plt.plot(np.abs(np.real(ideal_chirp_FFT)))
#  plt.title('Ideal Chirp Spectra - Amplitude (Standard Order)')
#  plt.subplot(4,1,4)
#  plt.plot(np.unwrap(np.angle(ideal_chirp_FFT)))
#  plt.title('Ideal Chirp Spectra - Unwrapped Phase (Standard Order)')
#  plt.xlabel('Sample')
#  plt.tight_layout()
#  plt.savefig('../runs/FritzFigsRequest/IdealChirp_FFT_Standard.png', dpi=1000)
  #
  # Load nal_filter.dat
  #
  cal_filter = np.fromfile('../calib/cal_filter.dat', '<f')
#  plt.clf()
#  plt.plot(cal_filter)
#  plt.title('Raw cal_filter')
#  plt.xlabel('Sample')
#  plt.savefig('../runs/FritzFigsRequest/cal_filter_raw.png', dpi=1000)
#  plt.show()
  cal_filter = cal_filter[:1800] + 1j*cal_filter[1800:]
 # plt.clf()
 # plt.subplot(4,1,1)
 # plt.plot(np.real(cal_filter))
 # plt.title('cal_filter - Real')
 # plt.subplot(4,1,2)
 # plt.plot(np.imag(cal_filter))
 # plt.title('cal_filter - Imaginary')
 # plt.subplot(4,1,3)
 # plt.plot(np.abs(np.real(cal_filter)))
 # plt.title('cal_filter - Amplitude')
 # plt.subplot(4,1,4)
 # plt.plot(np.unwrap(np.angle(cal_filter)))
 # plt.title('cal_filter - Unwrapped Phase')
 # plt.tight_layout()
 # plt.savefig('../runs/FritzFigsRequest/cal_filter_ComplexSpectra.png', dpi=1000)
  cal_filter = np.roll(cal_filter, 900)
#  plt.clf()
#  plt.subplot(4,1,1)
#  plt.plot(np.real(cal_filter))
#  plt.title('Rolled cal_filter - Real')
#  plt.subplot(4,1,2)
#  plt.plot(np.imag(cal_filter))
#  plt.title('Rolled cal_filter - Imaginary')
#  plt.subplot(4,1,3)
#  plt.plot(np.abs(np.real(cal_filter)))
#  plt.title('Rolled cal_filter - Amplitude')
#  plt.subplot(4,1,4)
#  plt.plot(np.unwrap(np.angle(cal_filter)))
#  plt.title('Rolled cal_filter - Unwrapped Phase')
#  plt.tight_layout()
#  plt.savefig('../runs/FritzFigsRequest/Rolled_cal_filter_ComplexSpectra.png', dpi=1000)
#  plt.show()
#  sys.exit()
  dechirp = np.zeros(3600, complex)
  dechirp[1800:] = cal_filter
  dechirp = (dechirp*ideal_chirp_FFT)
  dechirp_IFFT = np.fft.ifft(dechirp)
#  plt.clf()
#  plt.subplot(4,1,1)
#  plt.plot(np.real(dechirp))
#  plt.title('cal_chirp_FFT - Real')
#  plt.subplot(4,1,2)
#  plt.plot(np.imag(dechirp))
#  plt.title('cal_chirp_FFT - Imaginary')
#  plt.subplot(4,1,3)
#  plt.plot(np.abs(np.real(dechirp)))
#  plt.title('cal_chirp_FFT - Amplitude')
#  plt.subplot(4,1,4)
#  plt.plot(np.unwrap(np.angle(dechirp)))
#  plt.title('cal_chirp_FFT - Unwrapped Phase')
#  plt.tight_layout()
#  plt.savefig('../runs/FritzFigsRequest/cal_filter_FFT_ComplexSpectra.png', dpi=1000)
#  plt.show()
#  plt.clf()
#  plt.subplot(2,1,1)
#  plt.plot(np.real(dechirp_IFFT))
#  plt.title('IFFT cal_chirp_FFT - Real')
#  plt.subplot(2,1,2)
#  plt.plot(np.imag(dechirp_IFFT))
#  plt.title('IFFT cal_chirp_FFT - Imaginary')
#  plt.xlabel('Sample')
#  plt.tight_layout()
#  plt.savefig('../runs/FritzFigsRequest/cal_filter_FFT_IFFT_TimeSeries.png', dpi=1000)
#  plt.show() 
  #
  # Load a reference chirp example
  #
  raw_ref_chirp = np.fromfile('../calib/reference_chirp_p00tx_p00rx.dat', '<f')
  #ref_chirp = np.conj(ref_chirp[:2048] + 1j*ref_chirp[2048:])
  raw_ref_chirp = (raw_ref_chirp[:2048] + 1j*raw_ref_chirp[2048:])
  ref_chirp = np.zeros(4096, complex)
  ref_chirp[2048:] = raw_ref_chirp
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
  plt.show()
  #
  #
  #
  return

if __name__ == '__main__':
  main()
