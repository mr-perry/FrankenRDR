#
# Import other libraries
#
from os.path import isfile
from numpy import fromfile

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
  SSInstrMode = { 'SS01': {'Presum': 32, 'BitsPerSample': 8},
                'SS02': {'Presum': 28, 'BitsPerSample': 6},
                'SS03': {'Presum': 16, 'BitsPerSample': 4},
                'SS04': {'Presum': 8, 'BitsPerSample': 8},
                'SS05': {'Presum': 4, 'BitsPerSample': 6},
                'SS06': {'Presum': 2, 'BitsPerSample': 4},
                'SS07': {'Presum': 1, 'BitsPerSample': 8},
                'SS08': {'Presum': 32, 'BitsPerSample': 6},
                'SS09': {'Presum': 28, 'BitsPerSample': 4},
                'SS10': {'Presum': 16, 'BitsPerSample': 8},
                'SS11': {'Presum': 8, 'BitsPerSample': 6},
                'SS12': {'Presum': 4, 'BitsPerSample': 4},
                'SS13': {'Presum': 2, 'BitsPerSample': 8},
                'SS14': {'Presum': 1, 'BitsPerSample': 6},
                'SS15': {'Presum': 32, 'BitsPerSample': 4},
                'SS16': {'Presum': 28, 'BitsPerSample': 8},
                'SS17': {'Presum': 16, 'BitsPerSample': 6},
                'SS18': {'Presum': 8, 'BitsPerSample': 4},
                'SS19': {'Presum': 4, 'BitsPerSample': 8},
                'SS20': {'Presum': 2, 'BitsPerSample': 6},
                'SS21': {'Presum': 1, 'BitsPerSample': 4},
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
      EDRData['OpMode'] = SSInstrMode[EDR['OpMode']]
      print(EDRData['OpMode'])
  return
    
    
def loadData(_file):
  if isfile(_file):
    data = fromfile(_file, dtype='<f') 
   return 
