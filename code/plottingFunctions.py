from scipy.signal import welch
import matplotlib.pyplot as plt
import numpy as np

def bytescl(array, mindata=None, maxdata=None, top=255):
  #
  # Byte scaling algorithm for greater contrast
  #
  if mindata is None: mindata = np.nanmin(array)
  if maxdata is None: maxdata = np.nanmax(array)
  scl = np.maximum(np.minimum(((top+0.9999)*(array-mindata)/(maxdata-mindata)).astype(np.int16), top),0)
  return scl

def rdr2san(data, fname='rdr2san', maxdb=0, top=255):
  pow_out = np.power(np.abs(data), 2)			# Convert data to power
  db = 10 * np.log10(pow_out)				# decibels 
  maxdb = np.amax(db)
  sig = db/maxdb*255
  sig[np.where(sig < 0)] = 0.				# Zero out values below noise floor
  sig[np.where(sig > 255)] = 255.			# Clip values greater than MAXDB
  imName = '../runs/' + str(fname) + '.eps'
  plt.imshow(sig, cmap='gray')
  plt.savefig(imName, format='eps', dpi=1000)
  plt.show()
  return

def plotEDR(data, fname='plotEDR_new', ptype='Amp', thres=0.3, rel=False):
  if ptype == 'Amp':
    bmpName = '../runs/' + str(fname) + '_amp.bmp'
    pngName = '../runs/' + str(fname) + '_amp.png'
    #
    # I use the modulus amplitude here
    #
    Amp = np.abs(data)
    if rel == True:
      mx = np.amax(Amp, axis=0)
    else:
      mx = np.amax(Amp)
    pic = Amp / mx
    pic[np.where(pic < thres)] = 0.0
  elif ptype == 'Pow':
    bmpName = '../runs/' + str(fname) + '_pow.bmp'
    pngName = '../runs/' + str(fname) + '_pow.png'
    Pow = np.power(np.abs(data),2)
    if rel == True:
      mx = np.amax(Pow, axis=0)
    else:
      mx = np.amax(Pow)
    pic = Pow / mx
    pic[np.where(pic < thres)] = 0.0
  elif ptype == 'dB':
    if thres > 0:
      print('Threshold set to high, setting to -20 dB')
      thres = -15
    bmpName = '../runs/' + str(fname) + '_dB.bmp'
    pngName = '../runs/' + str(fname) + '_dB.png'
    Pow = np.power(np.abs(data),2)
    if rel == True:
      mx = np.amax(Pow, axis=0)
    else:
      mx = np.amax(Pow)
    pic = 10*np.log10(Pow/mx)
    pic[np.where(pic < thres)] = -99.0
  plt.imshow(pic, cmap='gray', vmin=thres)
  plt.savefig(pngName, dpi=1000)
  plt.imsave(bmpName, pic, cmap='gray') 
  return


def makePSD(data, oneside=False, rel=True, dt=37.5e-9, title='PSD', fname='PSD', dpi=500):
  dt = dt
  f, Pxx = welch(data, 1/dt, return_onesided=oneside)
  f = np.fft.fftshift(f)
  Pxx = np.fft.fftshift(Pxx)
  if rel:
    Pxx = Pxx / np.amax(Pxx)
  plt.plot(f/1e6, Pxx)
  plt.xlabel('Frequency (MHz)')
  plt.ylabel('Power Relative to Max')
  plt.title(title)
  plt.tight_layout()
  pngName = fname + '.png'
  plt.savefig(pngName, dpi=dpi)
  return


def plotFirstReturn(data, type='Amp', sidelobe=False, title='First Return', fname='FirstReturn', dpi=500):
  if type == 'Amp':
    data = np.abs(np.real(data))
    plt.ylabel('Amplitude')
  if type == 'Real':
    data = np.real(data)
  elif type == 'Pow':
    data = np.power(np.abs(np.real(data)),2)
  elif type == 'dB':
    data = 10 * np.log10(np.power(np.abs(np.real(data))/np.amax(np.abs(np.real(data)))),2)
  idx = np.where(data == np.amax(data))
  t0 = idx[0] * 0.0375
  t = np.arange(idx[0]-10, idx[0]+20) * 0.0375
  t = t - t0
  plt.clf()
  plt.plot(t,data[int(idx[0]-10):int(idx[0]+20)])
  if sidelobe:
    plt.axvline(0.24)
  plt.xlabel('Time (us)')
  plt.title(title)
  pngName = fname + '.png'
  plt.savefig(pngName, dpi=dpi)
  return
