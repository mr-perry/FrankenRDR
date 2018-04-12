from os.path import isfile
import sys
def loadLBL(fname):
  if isfile(fname):
    lblFile = open(fname, 'r')
  else:
    print('{} file not found. Please check path and try again.'.format(fname))
  return lblFile

def parseLBL(_file):
  #
  # This function parses a SHARAD EDR LBL file and returns
  # a Python dictionary object containing the information
  #
  # Define dictionary
  #
  lblDic = {  'INSTR_MODE_ID': [],
              'PRI': [],
           }
  lines = _file.readlines()
  #
  # Remove end of line characters from list
  #
  lines = [x.rstrip('\n') for x in lines]
  lineCount = len(lines)
  #
  # Remove all blank rows
  #
  lines = [x for x in lines if x]
  print("{} empty lines removed.".format(lineCount - len(lines)))
  lineCount = len(lines)
  #
  # Remove comments
  #
  lines = [x for x in lines if "/*" not in x]
  print("{} comment lines removed.".format(lineCount - len(lines)))
  lineCount = len(lines)
  #
  # Start parsing
  #
  print("Parsing {} lines in LBL file.".format(lineCount))
  for _i in range(lineCount):
    #
    # For this simple test all I should need out of the LBL file are:
    #  INSTRUMENT_MODE_ID
    #  
    if lines[_i].split('=')[0].strip() == 'INSTRUMENT_MODE_ID':
      lblDic['INSTR_MODE_ID'] = lines[_i].split('=')[1].strip() 
    if lines[_i].split('=')[0].strip() == 'MRO:PULSE_REPETITION_INTERVAL':
      lblDic['PRI'] = lines[_i].split('=')[1].strip() 
    if lines[_i].split('=')[0].strip() == 'MRO:MANUAL_GAIN_CONTROL':
      lblDic['GAIN_CONTROL'] = lines[_i].split('=')[1].strip() 
    if lines[_i].split('=')[0].strip() == 'MRO:COMPRESSION_SELECTION_FLAG':
      lblDic['COMPRESSION'] = lines[_i].split('=')[1].strip() 
  return lblDic 
