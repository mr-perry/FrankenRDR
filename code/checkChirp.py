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
  #
  # Load nal_filter.dat
  #
  cal_filter = np.fromfile('../calib/cal_filter.dat', '<f')
  cal_filter = cal_filter[:1800] + 1j*cal_filter[1800:]
  cal_filter = np.roll(cal_filter, 900)
  dechirp = np.zeros(3600, complex)
  dechirp[1800:] = cal_filter
  dechirp = np.conj(dechirp*ideal_chirp_FFT)
  #
  # Load a reference chirp example
  #
  ref_chirp = np.fromfile('../calib/reference_chirp_p00tx_p00rx.dat', '<f')
  #ref_chirp = np.conj(ref_chirp[:2048] + 1j*ref_chirp[2048:])
  ref_chirp = (ref_chirp[:2048] + 1j*ref_chirp[2048:])
  #
  # Let's plot
  #
  plt.subplot(3,2,1)
  plt.plot(np.real(np.fft.ifft(np.fft.ifftshift(dechirp))))
  plt.subplot(3,2,3)
  plt.plot(np.power(np.abs(dechirp[1800:]),2))
  plt.subplot(3,2,5)
  plt.plot(np.unwrap(np.angle(dechirp[1800:])))
  #
  # Let's plot
  #
  plt.subplot(3,2,2)
  plt.plot((np.real(np.fft.ifft(np.fft.ifftshift(ref_chirp)))))
  plt.subplot(3,2,4)
  plt.plot(np.power(np.abs(ref_chirp),2))
  plt.subplot(3,2,6)
  plt.plot(np.unwrap(np.angle(ref_chirp)))
  plt.show()
  #
  #
  #
  return

if __name__ == '__main__':
  main()
