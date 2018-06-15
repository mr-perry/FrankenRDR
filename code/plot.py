import numpy as np
from plottingFunctions import bytescl, plotEDR, rdr2san

def main():
  fname = '../runs/td7_beta0_ps4.npy'
  data = np.load(fname)
#  rdr2san(data, fname='td4_test_shift_Fritz_NW_nofinalpad')
  plotEDR(data, fname='td7_beta0_ps4', rel=True)
  return

if __name__ == '__main__':
  main()
