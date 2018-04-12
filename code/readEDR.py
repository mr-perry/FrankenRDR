#
# Import other libraries
#
import sys
from os.path import isfile, basename
from numpy import fromfile, dtype

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
    
    
def loadData(_file, OpMode):
  if isfile(_file):
    if OpMode['BitsPerSample'] == 4:
      dtype = '>i4'
    elif OpMode['BitsPerSample'] == 6:
      dtype = '>i6'
    elif OpMode['BitsPerSample'] == 8:
      dtype = '>i8'
    print(OpMode['BitsPerSample'])
    data = fromfile(_file, dtype=dtype) 
    print(len(data))
  else:
    print('File not found')
    sys.exit()
  return 

def main():
  fname = '../data/e_0168901_001_ss19_700_a_s.dat'
  EDRData = parseEDRname(fname)
  loadData(fname, EDRData['OpMode'])
  return

if __name__ == '__main__':
  main()
