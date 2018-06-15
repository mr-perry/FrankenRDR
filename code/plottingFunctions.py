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
#        print,'Minimum byte-scaled signal-above-noise is ',min(sig)
#        print,'Maximum byte-scaled signal-above-noise is ',max(sig)
  imName = '../runs/' + str(fname) + '.eps'
  plt.imshow(sig, cmap='gray')
  plt.savefig(imName, format='eps', dpi=1000)
  plt.show()
  return

def plotEDR(data, fname='plotEDR', rel=True, pmin=0.0, pmax=40):
  pmax = 0
  pmin = -20
  pow_out = np.power(np.abs(data), 2) 
  if rel:
    mx = np.max(pow_out, axis=0)
  else:
    mx = np.amax(pow_out)
  db = 10 * np.log10( pow_out / mx )					# this is how Bruce calcs. dB
#  f = db.flatten()
#  plt.hist(f, bins=100)
#  plt.yscale('log')
#  plt.show()
#  sys.exit()
  pic = bytescl(db, mindata = pmin, maxdata = pmax)  
  bmpName = '../runs/' + str(fname) + '.bmp'
  imName = '../runs/' + str(fname) + '.eps'
  plt.imshow(pic, cmap='gray')
  plt.imsave(bmpName, pic, cmap='gray')
#  plt.savefig(imName, format='eps', dpi=2000)
  plt.show()
  return


def plotEDR_new(data, fname='plotEDR_new'):
  bmpName = '../runs/' + str(fname) + '.bmp'
  Amp = np.abs(np.real(data))
  mx = np.amax(Amp, axis=0)
  mn = np.amin(Amp, axis=0)
  pic = np.zeros(np.shape(Amp))
  for _n in range(0, np.shape(Amp)[1]):
    trace = Amp[:, _n] / mx[_n]
    trace[np.where(trace < 0.3)] = 0.0
    pic[:, _n] = trace 
  plt.imshow(pic, cmap='gray')
  plt.imsave(bmpName, pic, cmap='gray') 
  plt.show()
  return

def main():
#  fname = '../runs/trueCentr.npy'
#  fname = '../runs/processEDR5.npy'
  fname = '../runs/processEDR4.npy'
  data = np.load(fname)
#  plotEDR(data, pmin=0, pmax=40)
#  rdr2san(data, maxdb=40)
  return

if __name__ == '__main__':
  main()
