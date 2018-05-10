import numpy as np
from plottingFunctions import bytescl, plotEDR, rdr2san

def main():
  fname = '../runs/WindowShift_KB6_32.npy'
  data = np.load(fname)
#  rdr2san(data, fname='NoWindowShift_KB6_rdr2san')
  plotEDR(data, fname='NoWindowShift_KB6', rel=False, pmin=-30, pmax=-20)
  return

if __name__ == '__main__':
  main()
