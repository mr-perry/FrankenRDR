#
# Import other libraries
#
import sys
from os.path import isfile, basename
from numpy import fromfile, dtype
import struct

def parseEDRname(fname):
  """
    This function parses the SHARAD EDR data file name and creates a
    dictionary of relevent information. The file names are organized by
    7 data fields all separated by '_'.
    Data Product: _file[0] (E for EDR)
    Transaction ID: _file[1:7] (Orbit number) 
                    _file[7:8] (number of the Operation Sequence Table (OST) that was used)
  """
  #
  # Set up dictionary for Instrument node
  #
  # Subsurface sounding
  #
  SSInstrMode = { 'SS01': {'Mode': 'SS01', 'Presum': 32, 'BitsPerSample': 8},
                'SS02': {'Mode': 'SS02','Presum': 28, 'BitsPerSample': 6},
                'SS03': {'Mode': 'SS03','Presum': 16, 'BitsPerSample': 4},
                'SS04': {'Mode': 'SS04','Presum': 8, 'BitsPerSample': 8},
                'SS05': {'Mode': 'SS05','Presum': 4, 'BitsPerSample': 6},
                'SS06': {'Mode': 'SS06','Presum': 2, 'BitsPerSample': 4},
                'SS07': {'Mode': 'SS07','Presum': 1, 'BitsPerSample': 8},
                'SS08': {'Mode': 'SS08','Presum': 32, 'BitsPerSample': 6},
                'SS09': {'Mode': 'SS09','Presum': 28, 'BitsPerSample': 4},
                'SS10': {'Mode': 'SS10','Presum': 16, 'BitsPerSample': 8},
                'SS11': {'Mode': 'SS11','Presum': 8, 'BitsPerSample': 6},
                'SS12': {'Mode': 'SS12','Presum': 4, 'BitsPerSample': 4},
                'SS13': {'Mode': 'SS13','Presum': 2, 'BitsPerSample': 8},
                'SS14': {'Mode': 'SS14','Presum': 1, 'BitsPerSample': 6},
                'SS15': {'Mode': 'SS15','Presum': 32, 'BitsPerSample': 4},
                'SS16': {'Mode': 'SS16','Presum': 28, 'BitsPerSample': 8},
                'SS17': {'Mode': 'SS17','Presum': 16, 'BitsPerSample': 6},
                'SS18': {'Mode': 'SS18','Presum': 8, 'BitsPerSample': 4},
                'SS19': {'Mode': 'SS19','Presum': 4, 'BitsPerSample': 8},
                'SS20': {'Mode': 'SS20','Presum': 2, 'BitsPerSample': 6},
                'SS21': {'Mode': 'SS21','Presum': 1, 'BitsPerSample': 4},
               }
  #
  # Receive only
  #
  ROInstrMode = { 'RO01': {'Presum': 32, 'BitsPerSample': 8},
                'RO02': {'Presum': 28, 'BitsPerSample': 6},
                'RO03': {'Presum': 16, 'BitsPerSample': 4},
                'RO04': {'Presum': 8, 'BitsPerSample': 8},
                'RO05': {'Presum': 4, 'BitsPerSample': 6},
                'RO06': {'Presum': 2, 'BitsPerSample': 4},
                'RO07': {'Presum': 1, 'BitsPerSample': 8},
                'RO08': {'Presum': 32, 'BitsPerSample': 6},
                'RO09': {'Presum': 28, 'BitsPerSample': 4},
                'RO10': {'Presum': 16, 'BitsPerSample': 8},
                'RO11': {'Presum': 8, 'BitsPerSample': 6},
                'RO12': {'Presum': 4, 'BitsPerSample': 4},
                'RO13': {'Presum': 2, 'BitsPerSample': 8},
                'RO14': {'Presum': 1, 'BitsPerSample': 6},
                'RO15': {'Presum': 32, 'BitsPerSample': 4},
                'RO16': {'Presum': 28, 'BitsPerSample': 8},
                'RO17': {'Presum': 16, 'BitsPerSample': 6},
                'RO18': {'Presum': 8, 'BitsPerSample': 4},
                'RO19': {'Presum': 4, 'BitsPerSample': 8},
                'RO20': {'Presum': 2, 'BitsPerSample': 6},
                'RO21': {'Presum': 1, 'BitsPerSample': 4},
               }
  if type(fname) is str:
    #
    # Get base file name
    #
    fname = basename(fname)
    fname = fname.split("_")
    EDRData = {'DataProduct': fname[0],
               'OrbitNum': fname[1][0:5],
               'OSTnum': fname[1][5:],
               'OSTLineNum': fname[2],
               'OpMode': fname[3].upper(),
               'PRF': fname[4],
               'Version': fname[5],
               'ftype': fname[6].split('.')[0],
               'ext': fname[6].split('.')[1]
              }
    if EDRData['OpMode'][0:2] == 'SS':
      EDRData['OpMode'] = SSInstrMode[EDRData['OpMode']]
    elif EDRData['OpMode'][0:2] == "RO":
      EDRData['OpMode'] = ROInstrMode[EDRData['OpMode']]
  return EDRData
    
    
def readBlock(_file, rowNum, BitsPerSample):
  """
    This function will read in a block of data from the science telemetry
    file and return the Ancilliary and Echo data for that block
  """
  #
  # Determine record length
  #
  if BitsPerSample == 4: 
    recLen = 1986
  if BitsPerSample == 6: 
    recLen = 2886
  if BitsPerSample == 8: 
    recLen = 3786
  #
  # See if file was a path or an object
  #
  if type(_file) is str:
    _file = open(_file, 'rb')
  #
  # Point to the appropriate row within the binary data and read the record
  #
  _file.seek(rowNum*recLen)
  rawData = _file.read(recLen)
  #
  # Split data; rawA - raw ancilliary
  #  rawS - raw science data
  #
  rawA = rawData[0:186] 
  rawS = rawData[186:]
  #
  # Set up Ancilliary Dictionary
  #
  AncilData = { 'FMT_LENGTH': struct.unpack(">h",  rawA[10:12])[0]}
  #
  # Decompress science data
  #
  return AncilData, rawS
  

def decompressData(rawScience, Compression, presum, bit_resolution):
  """
    This function decompresses the raw science data based on the
    method of compression used as described in the ancilliary data file
    as well as the PDS Label file
  """
  if Compression == 'Static' or Compression == 'Fixed':
    """
      For fixed scaling, the uncompressed value U of the echo sample is given by:
        U = C 2^s / N
      where
        C is the compressed value
        S is given by L - R + 8
        N is the number of presummed echoes
        L is the base 2 logarithm of N rounded to the nearest integer
        R is the bit resolution of the compressed value (i.e. 4,6,8)
    """
    L = int(np.log2(N)) 
    S = L - bit_resolution + 8
    decompress
def main():
  fname = '../data/e_0168901_001_ss19_700_a_s.dat'
  EDRData = parseEDRname(fname)
  for _i in range(100):
    readBlock(fname, _i, 8)
  return

if __name__ == '__main__':
  main()
