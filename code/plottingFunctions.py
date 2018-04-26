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

def rdr2san(data, maxdb=40, top=255):
  pow_out = np.power(np.abs(data), 2)			# Convert data to power
  db = 10 * np.log10(pow_out)				# decibels 
   
  sig = db/maxdb*255
  sig[np.where(sig < 0)] = 0.				# Zero out values below noise floor
  sig[np.where(sig > 255)] = 255.			# Clip values greater than MAXDB
#        print,'Minimum byte-scaled signal-above-noise is ',min(sig)
#        print,'Maximum byte-scaled signal-above-noise is ',max(sig)
  plt.imshow(sig, cmap='gray')
  plt.savefig('../runs/rdr2san.eps', format='eps', dpi=1000)
  return

def plotEDR(data, pmin=0.0, pmax=40):
  pow_out = np.power(np.abs(data), 2)
  db = 10 * np.log10( pow_out )					# this is how Bruce calcs. dB
  pic = bytescl(db, mindata = pmin, maxdata = pmax)  
  plt.imshow(pic, cmap='gray')
  plt.savefig('../runs/processEDR.eps', format='eps', dpi=1000)
  return

def main():
#  fname = '../runs/trueCentr.npy'
#  fname = '../runs/processEDR5.npy'
  fname = '../runs/processEDR4.npy'
  data = np.load(fname)
  plotEDR(data, pmin=0, pmax=40)
  rdr2san(data, maxdb=40)
  return

if __name__ == '__main__':
  main()
