import numpy as np
from SARProcessing import rangeCompression
from readChirp import detChirpFiles

def decompressionTable():
  look_i = np.zeros(256,np.float)
  for ii in range(len(look_i)):
    #
    # Zero to positive maximum value
    #
    look_i[ii] = np.float(ii) if ( ii <= 127) else np.float(ii) - 256.00
  return look_i


def main():
  file1 = '../data/DEC_DATA/OBS_323403000_1_Mode_001.dat'
  look_i = decompressionTable()
  decoder_data = np.fromfile(file1, dtype=np.uint8, count=-1)
  nrec = 5602
  recLen = 3600
  decoder_data = decoder_data.reshape([recLen, nrec])
  for _n in range(0, nrec):
    range_record = look_i[decoder_data[:,_n]]
  print(np.amin(range_record), np.amax(range_record))
  return

if __name__ == '__main__':
  main()
