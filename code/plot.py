import numpy as np
from plottingFunctions import bytescl, plotEDR, rdr2san

def main():
  fname = '../runs/processEDR4.npy'
  data = np.load(fname)
  rdr2san(data)
  return

if __name__ == '__main__':
  main()
