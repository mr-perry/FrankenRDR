#
# Import Libraries
#
from scipy import signal
import time
import pandas as pd
import numpy as np
import struct
import bitstring
#from readLBL import loadLBL, parseLBL
#from readAuxFile import readAuxFile, readRow, detCalibChirp, loadCalChirp
#from readEDR import parseEDRname, readBlock
import matplotlib.pyplot as plt
import glob, os, sys


def parseLBLFile(fname):
  #
  # Define instrument modes
  #
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
                  'RO21': {'Presum': 1, 'BitsPerSample': 4}
                }
  #
  # Now parse LBL File
  #
  if os.path.isfile(fname):
    #
    # Initialize dictionary
    #
    lblDic = {'INSTR_MODE_ID': [],
              'PRI': [],
              'GAIN_CONTROL': [],
              'COMPRESSION': [],
              'RECORD_BYTES': [],
              'FILE_RECORDS': []
              }
    with open(fname) as f:
      lines = f.readlines() 
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
      if lines[_i].split('=')[0].strip() == 'RECORD_BYTES':
        if lblDic['RECORD_BYTES'] == []:
          lblDic['RECORD_BYTES'] = int(lines[_i].split('=')[1].strip())
      if lines[_i].split('=')[0].strip() == 'FILE_RECORDS':
        if lblDic['FILE_RECORDS'] == []:
          lblDic['FILE_RECORDS'] = int(lines[_i].split('=')[1].strip())
    #
    # Find the instrument mode
    #
    if lblDic['INSTR_MODE_ID'][0:2] == 'SS':
      lblDic['INSTR_MODE_ID'] = SSInstrMode[lblDic['INSTR_MODE_ID']]
    elif lblDic['INSTR_MODE_ID'][0:2] == 'RO':
      lblDic['INSTR_MODE_ID'] = ROInstrMode[lblDic['INSTR_MODE_ID']]
    return lblDic
  else:
    print("{} file not found. Please check path and try again.".format(fname))
    return


def parseAuxFile(fname):
  """
  length from label filke
  This function will read in the auxilliary file and return a pandas
  dataframe with the necessary data
  dt = np.dtype([('SCET_BLOCK_WHOLE', 'u4'),
                ('SCET_BLOCK_FRAC', 'u2'),
                ('EPHEMERIS_TIME', 'f8'),
                ('GEOMETRY_EPOCH', 'U23'),
                ('SOLAR_LONGITUDE', 'f8'),
                ('ORBIT_NUMBER', 'i4'),
                ('X_MARS_SC_POSITION_VECTOR', 'f8'),
                ('Y_MARS_SC_POSITION_VECTOR', 'f8'),
                ('Z_MARS_SC_POSITION_VECTOR', 'f8'),
                ('SPACECRAFT_ALTITUDE', 'f8'),
                ('SUB_SC_EAST_LONGITUDE', 'f8'),
                ('SUB_SC_PLANETOCENTRIC_LATITUDE', 'f8'),
                ('SUB_SC_PLANETOGRAPHIC_LATITUDE', 'f8'),
                ('X_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('Y_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('Z_MARS_SC_VELOCITY_VECTOR', 'f8'),
                ('MARS_SC_RADIAL_VELOCITY', 'f8'),
                ('MARS_SC_TANGENTIAL_VELOCITY', 'f8'),
                ('LOCAL_TRUE_SOLAR_TIME', 'f8'),
                ('SOLAR_ZENITH_ANGLE', 'f8'),
                ('SC_PITCH_ANGLE', 'f8'),
                ('SC_YAW_ANGLE', 'f8'),
                ('SC_ROLL_ANGLE', 'f8'),
                ('MRO_SAMX_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAMX_OUTER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAPX_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_SAPX_OUTER_GIMBAL_ANGLE', 'f8'),
                ('MRO_HGA_INNER_GIMBAL_ANGLE', 'f8'),
                ('MRO_HGA_OUTER_GIMBAL_ANGLE', 'f8'),
                ('DES_TEMP', 'f4'),
                ('DES_5V', 'f4'),
                ('DES_12V', 'f4'),
                ('DES_2V5', 'f4'),
                ('RX_TEMP', 'f4'),
                ('TX_TEMP', 'f4'),
                ('TX_LEV', 'f4'),
                ('TX_CURR', 'f4'),
                ('CORRUPTED_DATA_FLAG', 'i2')
               ]
               )
  """
  #
  # Set up dictionary
  #
  auxData ={'SCET_BLOCK_WHOLE': [],
            'SCET_BLOCK_FRAC': [],
            'EPHEMERIS_TIME': [],
            'ELAPSED_TIME': [],
            'GEOMETRY_EPOCH': [],
            'SOLAR_LONGITUDE': [],
            'ORBIT_NUMBER': [],
            'X_MARS_SC_POSITION_VECTOR': [],
            'Y_MARS_SC_POSITION_VECTOR': [],
            'Z_MARS_SC_POSITION_VECTOR': [],
            'SPACECRAFT_ALTITUDE': [],
            'SUB_SC_EAST_LONGITUDE': [],
            'SUB_SC_PLANETOCENTRIC_LATITUDE': [],
            'SUB_SC_PLANETOGRAPHIC_LATITUDE': [],
            'X_MARS_SC_VELOCITY_VECTOR': [],
            'Y_MARS_SC_VELOCITY_VECTOR': [],
            'Z_MARS_SC_VELOCITY_VECTOR': [],
            'MARS_SC_RADIAL_VELOCITY': [],
            'MARS_SC_TANGENTIAL_VELOCITY': [],
            'LOCAL_TRUE_SOLAR_TIME': [],
            'SOLAR_ZENITH_ANGLE': [],
            'SC_PITCH_ANGLE': [],
            'SC_YAW_ANGLE': [],
            'SC_ROLL_ANGLE': [],
            'MRO_SAMX_INNER_GIMBAL_ANGLE': [],
            'MRO_SAMX_OUTER_GIMBAL_ANGLE': [],
            'MRO_SAPX_INNER_GIMBAL_ANGLE': [],
            'MRO_SAPX_OUTER_GIMBAL_ANGLE': [],
            'MRO_HGA_INNER_GIMBAL_ANGLE': [],
            'MRO_HGA_OUTER_GIMBAL_ANGLE': [],
            'DES_TEMP': [],
            'DES_5V': [],
            'DES_12V': [],
            'DES_2V5': [],
            'RX_TEMP': [],
            'TX_TEMP': [],
            'TX_LEV': [],
            'TX_CURR': [],
            'CORRUPTED_DATA_FLAG': []
          }
  #
  # Each record is composed of 267 bytes
  #
  recLen = 267
  if os.path.isfile(fname):
    _file = open(fname, 'rb')
    fsize = os.path.getsize(fname)
    for _i in range(int(fsize/recLen)): # Go through all the rows
      _file.seek(_i*recLen)
      rawData = _file.read(recLen)
      auxData['SCET_BLOCK_WHOLE'].append(struct.unpack(">I", rawData[0:4])[0])
      auxData['SCET_BLOCK_FRAC'].append(struct.unpack(">H", rawData[4:6])[0])
      auxData['EPHEMERIS_TIME'].append(struct.unpack(">d", rawData[6:14])[0])
      auxData['ELAPSED_TIME'].append(auxData['EPHEMERIS_TIME'][_i] - auxData['EPHEMERIS_TIME'][0])
      auxData['GEOMETRY_EPOCH'].append(rawData[14:37].decode("utf-8"))
      auxData['SOLAR_LONGITUDE'].append(struct.unpack(">d", rawData[37:45])[0])
      auxData['ORBIT_NUMBER'].append(struct.unpack(">i", rawData[45:49])[0])
      auxData['X_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[49:57])[0])
      auxData['Y_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[57:65])[0])
      auxData['Z_MARS_SC_POSITION_VECTOR'].append(struct.unpack(">d", rawData[65:73])[0])
      auxData['SPACECRAFT_ALTITUDE'].append(struct.unpack(">d", rawData[73:81])[0])
      auxData['SUB_SC_EAST_LONGITUDE'].append(struct.unpack(">d", rawData[81:89])[0])
      auxData['SUB_SC_PLANETOCENTRIC_LATITUDE'].append(struct.unpack(">d", rawData[89:97])[0])
      auxData['SUB_SC_PLANETOGRAPHIC_LATITUDE'].append(struct.unpack(">d", rawData[97:105])[0])
      auxData['X_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[105:113])[0])
      auxData['Y_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[113:121])[0])
      auxData['Z_MARS_SC_VELOCITY_VECTOR'].append(struct.unpack(">d", rawData[121:129])[0])
      auxData['MARS_SC_RADIAL_VELOCITY'].append(struct.unpack(">d", rawData[129:137])[0])
      auxData['MARS_SC_TANGENTIAL_VELOCITY'].append(struct.unpack(">d", rawData[137:145])[0])
      auxData['LOCAL_TRUE_SOLAR_TIME'].append(struct.unpack(">d", rawData[145:153])[0])
      auxData['SOLAR_ZENITH_ANGLE'].append(struct.unpack(">d", rawData[153:161])[0])
      auxData['SC_PITCH_ANGLE'].append(struct.unpack(">d", rawData[161:169])[0])
      auxData['SC_YAW_ANGLE'].append(struct.unpack(">d", rawData[169:177])[0])
      auxData['SC_ROLL_ANGLE'].append(struct.unpack(">d", rawData[177:185])[0])
      auxData['MRO_SAMX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[185:193])[0])
      auxData['MRO_SAMX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[193:201])[0])
      auxData['MRO_SAPX_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[201:209])[0])
      auxData['MRO_SAPX_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[209:217])[0])
      auxData['MRO_HGA_INNER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[217:225])[0])
      auxData['MRO_HGA_OUTER_GIMBAL_ANGLE'].append(struct.unpack(">d", rawData[225:233])[0])
      auxData['DES_TEMP'].append(struct.unpack(">f", rawData[233:237])[0])
      auxData['DES_5V'].append(struct.unpack(">f", rawData[237:241])[0])
      auxData['DES_12V'].append(struct.unpack(">f", rawData[241:245])[0])
      auxData['DES_2V5'].append(struct.unpack(">f", rawData[245:249])[0])
      auxData['RX_TEMP'].append(struct.unpack(">f", rawData[249:253])[0])
      auxData['TX_TEMP'].append(struct.unpack(">f", rawData[253:257])[0])
      auxData['TX_LEV'].append(struct.unpack(">f", rawData[257:261])[0])
      auxData['TX_CURR'].append(struct.unpack(">f", rawData[261:265])[0])
      auxData['CORRUPTED_DATA_FLAG'].append(struct.unpack(">h", rawData[265:267])[0])
    #
    # Turn into Pandas dataframe because dictionaries are lame-sauce
    #
    df = pd.DataFrame.from_dict(auxData)
    df.to_csv('test.csv')
    #
    # Print out a summary of the data parsed
    #
    print(df.describe())
    #
    # Produce plots from Aux Data
    #
    # Use plt.plot rather than pandas wrapper
    #
#    makeAuxPlots(df)
    return df

def makeAuxPlots(df):
  f, axarr = plt.subplots(2, sharex=True)
  f.suptitle('Sharing X axis')
  X = df['ELAPSED_TIME']
  axarr[0].plot(X, df['SOLAR_LONGITUDE'], 'k.')
  axarr[1].plot(X, df['X_MARS_SC_POSITION_VECTOR'], 'k.') 
  axarr[1].plot(X, df['Y_MARS_SC_POSITION_VECTOR'], 'r.') 
  axarr[1].plot(X, df['Z_MARS_SC_POSITION_VECTOR'], 'b.') 
  plt.show()
  return


def parseAncilliary(data):
  """
    This function parses the ancilliary data for an individual SHARAD EDR data record
  """
  #
  # Check that length of data is 186 bytes
  #
  if len(data) != 186:
    print('Incorrect data supplied. Please ensure data is a 186 byte string')
    return
  else:
    #
    # Set up dictionary for items
    #
    #ancilliaryData = { 'SCET_BLOCK_WHOLE': struct.unpack('>I', data[0:4])[0],
    ancilliaryData = { 'SCET_BLOCK_WHOLE': bitstring.BitArray(data[0:4]).uint,
                       'SCET_BLOCK_FRAC': bitstring.BitArray(data[4:6]).uint,
                       'TLM_COUNTER': struct.unpack('>I', data[6:10])[0],
                       'FMT_LENGTH': struct.unpack('>H', data[10:12])[0],
                       'SPARE1': struct.unpack('>H', data[12:14])[0],
                       'SCET_OST_WHOLE': struct.unpack('>I', data[14:18])[0],
                       'SCET_OST_FRAC': struct.unpack('>H', data[18:20])[0],
                       'SPARE2': struct.unpack('>B', data[20:21])[0],
                       'OST_LINE_NUMBER': struct.unpack('>B', data[21:22])[0],
                       'OST_LINE': { },
                       'SPARE3': struct.unpack('>B', data[38:39])[0],
                       'DATA_BLOCK_ID': bitstring.BitArray(data[39:42]).uint,
                       'SCIENCE_DATA_SOURCE_COUNTER': struct.unpack('>H', data[42:44])[0],
                       'PACKET_SEGMENTATION_AND_FPGA_STATUS': { },
                       'SPARE4': struct.unpack('>B', data[46:47])[0], 
                       'DATA_BLOCK_FIRST_PRI': bitstring.BitArray(data[47:50]).uint,
                       'TIME_DATA_BLOCK_WHOLE': struct.unpack('>I', data[50:54])[0],
                       'TIME_DATA_BLOCK_FRAC': struct.unpack('>H', data[54:56])[0],
                       'SDI_BIT_FIELD': struct.unpack('>H', data[56:58])[0],
                       'TIME_N': struct.unpack('>f', data[58:62])[0],
                       'RADIUS_N': struct.unpack('>f', data[62:66])[0],
                       'TANGENTIAL_VELOCITY_N': struct.unpack('>f', data[66:70])[0],
                       'RADIAL_VELOCITY_N': struct.unpack('>f', data[70:74])[0],
                       'TLP': struct.unpack('>f', data[74:78])[0],
                       'TIME_WPF': struct.unpack('>f', data[78:82])[0],
                       'DELTA_TIME': struct.unpack('>f', data[82:86])[0],
                       'TLP_INTERPOLATE': struct.unpack('>f', data[86:90])[0],
                       'RADIUS_INTERPOLATE': struct.unpack('>f', data[90:94])[0],
                       'TANGENTIAL_VELOCITY_INTERPOLATE': struct.unpack('>f', data[94:98])[0],
                       'RADIAL_VELOCITY_INTERPOLATE': struct.unpack('>f', data[98:102])[0],
                       'END_TLP': struct.unpack('>f', data[102:106])[0],
                       'S_COEFFS': np.array([struct.unpack('>f', data[106:110]),
                                         struct.unpack('>f', data[110:114]),
                                         struct.unpack('>f', data[114:118]),
                                         struct.unpack('>f', data[118:122]),
                                         struct.unpack('>f', data[122:126]),
                                         struct.unpack('>f', data[126:130]),
                                         struct.unpack('>f', data[130:134]),
                                         struct.unpack('>f', data[134:138])
                                        ]),
                       'C_COEFFS': np.array([struct.unpack('>f', data[138:142]),
                                         struct.unpack('>f', data[142:146]),
                                         struct.unpack('>f', data[146:150]),
                                         struct.unpack('>f', data[150:154]),
                                         struct.unpack('>f', data[154:158]),
                                         struct.unpack('>f', data[158:162]),
                                         struct.unpack('>f', data[162:166])
                                        ]),
                       'SLOPE': struct.unpack('>f', data[166:170])[0],
                       'TOPOGRAPHY': struct.unpack('>f', data[170:174])[0],
                       'PHASE_COMPENSATION_STEP': struct.unpack('>f', data[174:178])[0],
                       'RECEIVE_WINDOW_OPENING_TIME': struct.unpack('>f', data[178:182])[0],
                       'RECEIVE_WINDOW_POSITION': struct.unpack('>f', data[182:186])[0],
                     }
                      
    #####################################################################################
    #
    # PACKET_SEGMENTATION_AND_FPGA_STATUS bit string
    #
    PSAFS = bitstring.BitArray(data[44:46])
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SCIENTIFIC_DATA_TYPE'] = PSAFS[0:1].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SEGMENTATION_FLAG'] = PSAFS[1:3].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE1'] = PSAFS[3:8].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['SPARE2'] = PSAFS[8:12].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['DMA_ERROR'] = PSAFS[12:13].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TC_OVERRUN'] = PSAFS[13:14].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['FIFO_FULL'] = PSAFS[14:15].uint
    ancilliaryData['PACKET_SEGMENTATION_AND_FPGA_STATUS']['TEST'] = PSAFS[15:16].uint
    #####################################################################################
    #
    # OST_LINE_NUMBER BIT STRING
    #
    OST_LINE = bitstring.BitArray(data[22:39])
    ancilliaryData['OST_LINE']['PULSE_REPETITION_INTERVAL'] = OST_LINE[0:4].uint
    ancilliaryData['OST_LINE']['PHASE_COMPENSATION_TYPE'] = OST_LINE[4:8].uint
    ancilliaryData['OST_LINE']['SPARE1'] = OST_LINE[8:10].uint
    ancilliaryData['OST_LINE']['DATA_LENGTH_TAKEN'] = OST_LINE[10:32].uint
    ancilliaryData['OST_LINE']['OPERATIVE_MODE'] = OST_LINE[32:40].uint
    ancilliaryData['OST_LINE']['MANUAL_GAIN_CONTROL'] = OST_LINE[40:48].uint
    ancilliaryData['OST_LINE']['COMPRESSION_SELECTION'] = OST_LINE[48:49].bool
    ancilliaryData['OST_LINE']['CLOSED_LOOP_TRACKING'] = OST_LINE[49:50].bool
    ancilliaryData['OST_LINE']['TRACKING_DATA_STORAGE'] = OST_LINE[50:51].bool
    ancilliaryData['OST_LINE']['TRACKING_PRE_SUMMING'] = OST_LINE[51:54].uint
    ancilliaryData['OST_LINE']['TRACKING_LOGIC_SELECTION'] = OST_LINE[54:55].uint
    ancilliaryData['OST_LINE']['THRESHOLD_LOGIC_SELECTION'] = OST_LINE[55:56].uint
    ancilliaryData['OST_LINE']['SAMPLE_NUMBER'] = OST_LINE[56:60].uint
    ancilliaryData['OST_LINE']['SPARE2'] = OST_LINE[60:61].uint
    ancilliaryData['OST_LINE']['ALPHA_BETA'] = OST_LINE[61:63].uint
    ancilliaryData['OST_LINE']['REFERENCE_BIT'] = OST_LINE[63:64].uint
    ancilliaryData['OST_LINE']['THRESHOLD'] = OST_LINE[64:72].uint
    ancilliaryData['OST_LINE']['THRESHOLD_INCREMENT'] = OST_LINE[72:80].uint
    ancilliaryData['OST_LINE']['SPARE3'] = OST_LINE[80:84].uint
    ancilliaryData['OST_LINE']['INITIAL_ECHO_VALUE'] = OST_LINE[84:87].uint
    ancilliaryData['OST_LINE']['EXPECTED_ECHO_SHIFT'] = OST_LINE[87:90].uint
    ancilliaryData['OST_LINE']['WINDOW_LEFT_SHIFT'] = OST_LINE[90:93].uint
    ancilliaryData['OST_LINE']['WINDOW_RIGHT_SHIFT'] = OST_LINE[93:96].uint
    ancilliaryData['OST_LINE']['SPARE4'] = OST_LINE[96:128].uint
  return ancilliaryData


def readEDRrecord(_file, record, recLen, bps):
  """
    This function reads a single EDR record and returns the ancilliary 
    and science data. Data are read and decoded based off their bit resolution
    Input:
      _file: File object identifier pointing to the data file containing the binary data
    Output:
      ancil: Data dictionary containing all the information parsed from the ancilliary data for the record
      echoes: Decoded science data
  """
  #
  # Make sure we are at the beginning of a record
  #
  _file.seek(0)
  #
  # Now fast forward to the appropriate location
  _file.seek(record*recLen)
  #
  # Step 1: Read in binary data
  #
  rawData = _file.read(recLen)
  #
  # Separate ancil and echo data
  #
  # Step 2: Read and parse ancilliary data
  #
  ancil = rawData[:186]
  #
  # Parse Ancilliary Data
  #
  ancDic = parseAncilliary(ancil)
  #
  # Step 3: Separate Actual echoes
  #
  echoes = rawData[186:]
  #
  # Okay, now decode the data
  # For 8-bit samples, there is a sample for every byte
  # For 6-bit samples, there are 2 samples for every 3 bytes
  # For 4-bit samples, there are 2 samples for every byte
  #
  # Making the vector have 4096 rows is to allow for proper decovolution of the
  # calibrated chirp
  #
  decoded_data = np.zeros(4096, int)
  #
  # Step 4: Decode the science data based on the bit resolution
  #
  # Break Byte stream into Bit Stream
  # This isn't necessary for the 8-bit resolution, but to keep the code
  # clean, split the bytes up regardless of resolution
  #
  b = bitstring.BitArray(echoes)
  for _j in range(0, len(b), BitsPerSample):
    decoded_data[int(_j/BitsPerSample)] = b[_j:_j+BitsPerSample].int 
  return decoded_data, ancDic

def decompressSciData(data, compression):
  """
    This function decompresses the data based off page 8 of the 
    SHALLOW RADAR EXPERIMENT DATA RECORD SOFTWARE INTERFACE SPECIFICATION.
    If the compression type is 'Static' (i.e. False):
       U = C*2**S / N
         C is the Compressed Data
         S = L - R + 8
            L is base 2 log of N rounded UP to the nearest integer
            R is the bit resolution of the compressed values
         N is the number of pre-summed echoes
    If the compression type is 'dynamic' (i.e. True):
      NOT WORKING YET
      U = C*2**S/N
        C is the compressed data
        N is the number of pre-summed echoes
        S = SDI for SDI <= 5
        S = SDI-6 for 5 < SDI <= 16
        S = SDI-16 for SDI > 16
          where SDI is the SDI_BIT_FIELD parameter from ANCILLIARY DATA
  """
  # Step 5: Decompress the data
  # Note: Only Static decompression works at this time
  # 
  if compression == False: # Static scaling
    L = np.ceil(np.log2(int(instrPresum)))
    R = BitsPerSample
    S = L - R + 8
    N = instrPresum
    decompressed_data = data * (np.power(2, S)/N)
    return decompressed_data
  elif compression == True:#dynamic scaling
    print('This is not yet available.')
    sys.exit()

def rangeCompression(fname, filtering='Match'):
  """
     The following steps must occur during range compression:
      1) A single record is read in
      2) The ancilliary data is parsed and decoded
      3) Parse the science data (echoes)
      4) Decode the echoes based on their bit resolution
      5) Decompress the echoes using either the static or dynamic methods described on page 
         8 of the SHALLOW RADAR EXPERIMENT DATA RECORD SOFTWARE INTERFACE SPECIFICATION
      6) Once record is decompressed determine which calibrated chirp file is needed
      7) Extend the length of a data record from 3600 to 4096 using zero padding
      8) Multiply this raw data by the complex exponential of 2*pi*i*F_c*t,
         where F_c is equal to (80/3 - 20) MHz, and t is a vector 4096 samples long containing
         linearly increasing values starting from zero, and each spaced 3/80 microseconds apart.
      9) Compute the FFT of the result
     10) Select only the central 2048 samples of the complex spectrum, that is from sample 1025
         to sample 3072 if spectrum are ordered for increasing frequency.
     11) Multiply the result by the complex conjugate of the reference chirp
     12) Compute the IFFT of the result
  """
  #
  # Set up some variables needed later
  #
  F_c = (80/3 - 20)  
  t = np.linspace(0, 4096*0.0375, num=4096) 
  compExp = 2*np.pi*1j*F_c*t
  #
  # Determine record length based off bits per sample
  #
  if BitsPerSample == 4:
    recLen = 1986
  elif BitsPerSample == 6:
    recLen = 2886
  elif BitsPerSample == 8:
    recLen = 3786
  #
  # Check if file exists and begin processing
  #
  if not os.path.isfile(fname):
    print('File does not exist...exiting')
  else:
    _file = open(fname, 'rb')
    fsize = os.path.getsize(fname)
    length = fsize/recLen
    #
    # Make empty array for EDR data
    #
    EDRData = np.zeros([2048, int(length)], complex)
    #
    # Begin range processing
    #
    for _i in range(int(length)): # Go through all the records
      #
      # Read in single record
      #
      sci_data, ancDic = readEDRrecord(_file, _i, recLen, BitsPerSample)
      #
      # Decompress the science data
      #
      echoes = decompressSciData(sci_data, ancDic['OST_LINE']['COMPRESSION_SELECTION'])
      #
      # Determine calibrated chirp
      #
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_i], AuxDF['RX_TEMP'][_i])
      #
      #  8) Multiply this raw data by the complex exponential of 2*pi*i*F_c*t,
      #     where F_c is equal to (80/3 - 20) MHz, and t is a vector 4096 samples long containing
      #     linearly increasing values starting from zero, and each spaced 3/80 microseconds apart.
      #
      decom_data_c = echoes * compExp
      #
      # Compute FFT
      #
      # The values in the result follow so-called “standard” order: If A = fft(a, n), 
      # then A[0] contains the zero-frequency term (the sum of the signal), which is 
      # always purely real for real inputs. Then A[1:n/2] contains the positive-frequency 
      # terms, and A[n/2+1:] contains the negative-frequency terms, in order of decreasingly negative frequency.
      #
      echoSpec = np.fft.fft(decom_data_c)
      #
      # Now reorder FFT so we can take out the central frequencies
      # 
      echoSpec = np.fft.fftshift(echoSpec)
      #
      # Select on the central 2048 frequencies
      # Does this need to be a more robust windowing functions???
      #
      echoSpec = echoSpec[1024:3072]
      #
      # Remove the calibrated chirp from the science data (i.e. Range compression)
      #
      if filtering == 'Match':
        #
        # Calculate the complex conjugate of the calibrated chirp
        #
        dechirp = np.conj(calChirp)
        #
        # Emply match filtering techn
        #
        dechirped = echoSpec * dechirp
      elif filtering == 'Inverse':
        dechirpEchoSpec = echoSpec / calChirp
      #
      # Now, window the resulting spectrum for sidelobe supression
      #
      #
      # Now, reorder the dechirped Spectrum and IFFT
      #
      EDRData[:, _i] = np.fft.ifft(np.fft.ifftshift(dechirped))
  return EDRData


def detChirpFiles(TxTemp, RxTemp):
  """
  This function determines the appropriate calibrated chirp file to use
  for range compression and returns the decoded calibrated chirp. 
  This is solely based off the temperatures of the TxTemp and RxTemp.
  """
  calibRoot = '../calib/'
  calibName = 'reference_chirp'
  ext = '.dat'
  TxCalNames = ['m20tx', 'm15tx', 'm10tx', 'm05tx',
                'p00tx', 'p20tx', 'p40tx', 'p60tx']
  RxCalNames = ['m20rx', 'p00rx', 'p20rx', 'p40rx',
                'p60rx']
  #
  # Define vectors for Tx and Rx temps
  #
  Tx = [-20, -15, -10, -5, 0, 20, 40, 60]
  Rx = [-20, 0, 20, 40, 60]
  calibChirpFiles = []
  TxDiff = []
  RxDiff = []
  #
  # Find distance
  #
  TxDiff[:] = [abs(x - TxTemp) for x in Tx]
  RxDiff[:] = [abs(x - RxTemp) for x in Rx]
  #
  # Find the indices of the closest Tx and Rx value
  #
  calibTx = TxCalNames[TxDiff.index(min(TxDiff))]
  calibRx = RxCalNames[RxDiff.index(min(RxDiff))]
  #
  # Construct File name
  #
  calChirpFile = calibRoot + calibName + '_' + \
                 TxCalNames[TxDiff.index(min(TxDiff))] + '_' + \
                 RxCalNames[RxDiff.index(min(RxDiff))] + ext
  if os.path.isfile(calChirpFile):
    calChirp = np.fromfile(calChirpFile, dtype='<f')
    calChirp = calChirp[:2048] + 1j*calChirp[2048:]
  else:
    print('Calibrated chirp file not found...exiting.')
    sys.exit()
  return calChirp

def main():
  #
  # Ground processing steps per SHALLOW RADAR REDUCED DATA RECORD SOFTWARE INTERFACE SPECIFICATION
  # document version 1.0 20 July 2007
  #
  # Decompression and Pre-summing
  # Range Processing
  # Azimuth Processing
  # Calibration
  # Ionospheric Correction
  # Time Alignment of Echoes
  #
  # Set up global variables
  # These will be items read from the LBL file that are necessary for proper
  # processing and decoding
  #
  global nrec, instrMode, presum_onboard, BitsPerSample, AuxDF 
  #
  # 
  # TEST DATA
  #
  td = 4
  fil = 'Match'
  fft_length = 1536
  compression = False;
  if td == 0:
    auxname = '../data/e_0168901_001_ss19_700_a_a.dat'
    lblname = '../data/e_0168901_001_ss19_700_a.lbl' 
    edrname = '../data/e_0168901_001_ss19_700_a_s.dat' 
  elif td == 1:
    auxname = '../data/e_1601301_001_ss19_700_a_a.dat'
    lblname = '../data/e_1601301_001_ss19_700_a.lbl' 
    edrname = '../data/e_1601301_001_ss19_700_a_s.dat' 
  elif td == 2:
    auxname = '../data/e_1408801_001_ss11_700_a_a.dat'
    lblname = '../data/e_1408801_001_ss11_700_a.lbl' 
    edrname = '../data/e_1408801_001_ss11_700_a_s.dat' 
  elif td == 3:
    auxname = '../data/e_0169401_001_ss18_700_a_a.dat'
    lblname = '../data/e_0169401_001_ss18_700_a.lbl' 
    edrname = '../data/e_0169401_001_ss18_700_a_s.dat' 
  elif td == 4: #VERY FLAT REGION
    auxname = '../data/e_0323403_001_ss19_700_a_a.dat'
    lblname = '../data/e_0323403_001_ss19_700_a.lbl' 
    edrname = '../data/e_0323403_001_ss19_700_a_s.dat' 
  #
  # Parse important data from label file
  #
  print('Reading label file...')
  lblDic = parseLBLFile(lblname)
  print('Finished reading label file...')
  nrec = lblDic['FILE_RECORDS']			# Number of records
  instrMode = lblDic['INSTR_MODE_ID']['Mode']
  presum_onboard = lblDic['INSTR_MODE_ID']['Presum']
  BitsPerSample = lblDic['INSTR_MODE_ID']['BitsPerSample']
  #
  # Okay, now that the label file has been parsed for the necessary information
  # read in the entire Aux file
  print('Reading Auxilliary file...')
  AuxDF = parseAuxFile(auxname)  
  print('Finished reading Auxilliary file...')
  #############################################################################
  #
  #
  n_range = 3600					# Typical SHARAD record length
  timing_error = 12.5 / 1.0e6				#12.5-us timing error; I don't know if I'll actually need this
  delay_res = 135.0e-06 / 3600.				# Delay resolution (sec)
  sharad_ipp = 1. / 700.28				# Pulse interval (sec)
  nflat = 49152 / presum_onboard			# No. of raw pulses at presum-1
  presum_proc = 16 if fft_length <= 1536 else 8		# Set internal presumming to balance the speed with radargram
 							# aliasing. Really should not push the aperture length much
							# past 3072 in presum-4 units (18-seconds). Presum factor 
							# 16X for PDS default, 8X for "long" mode.
  nsums = np.copy(nrec)
  #
  # Determine record length based off bits per sample
  #
  if BitsPerSample == 4:
    recLen = 1986
  elif BitsPerSample == 6:
    recLen = 2886
  elif BitsPerSample == 8:
    recLen = 3786
  ############################################################################
  #
  # At this point in the FPB, Bruce reads the header file and constructs
  # Vx and Vz polynomials...do I need to do this?
  #
  # Then he reads maptrack.out for information, can I get this from ancilliary data? 
  #
  # 
  #
  # Now he makes a windowing function, but zeros out half of the spectrum...I need to look through
  # some literature about this
  #
  ############################################################################
  print('Reading, decoding, and decompressing science telemetry file...')
  if nflat >= nsums:
    nflat = nsums - 1 
  fft_length_pass = nflat
  #
  # Check if EDR file exists
  #
  if not os.path.isfile(edrname):
    print('File does not exist...exiting')
  else:
    _file = open(edrname, 'rb')
  #
  # All of the input pulses in groups
  #
  nsample = -1
  for _k in range(0, nsums, int(presum_proc/presum_onboard)):
    nsample += 1
    if _k+(presum_proc/presum_onboard-1) < nsums:
      presum_record = np.zeros(2048, complex)		# Zero out this stack
      range_record = np.zeros(4096, float)		# Zero out this stack
      #
      # Point to the beginning of this data block in the raw data file
      #
      #
      # Make sure we are at the beginning of a record
      #
      _file.seek(0)
      #
      # Now fast forward to the appropriate location
      _file.seek(_k*recLen)
      for _m in range(0, int(presum_proc/presum_onboard)):
        rawData = _file.read(recLen)
        #ancil = rawData[:186]
        #
        # Parse Ancilliary Data
        #
        #ancDic = parseAncilliary(ancil)
        #
        # Step 3: Separate Actual echoes
        #
        echoes = rawData[186:]
        #
        # Decode the echoes
        #
        #
        # Okay, now decode the data
        # For 8-bit samples, there is a sample for every byte
        # For 6-bit samples, there are 2 samples for every 3 bytes
        # For 4-bit samples, there are 2 samples for every byte
        #
        # Making the vector have 4096 rows is to allow for proper decovolution of the
        # calibrated chirp
        #
        decoded_data = np.zeros(4096, int)
        #
        # Step 4: Decode the science data based on the bit resolution
        #
        # Break Byte stream into Bit Stream
        # This isn't necessary for the 8-bit resolution, but to keep the code
        # clean, split the bytes up regardless of resolution
        #
        b = bitstring.BitArray(echoes)
        for _j in range(0, len(b), BitsPerSample):
          decoded_data[int(_j/BitsPerSample)] = b[_j:_j+BitsPerSample].int 
        #
        # Now decompress the data
        #
        if compression == False: # Static scaling
          L = np.ceil(np.log2(int(presum_onboard)))
          R = BitsPerSample
          S = L - R + 8
          N = presum_onboard
          decompressed_data = decoded_data * (np.power(2, S)/N)
        elif compression == True:#dynamic scaling
          print('This is not yet available.')
          sys.exit()
        range_record = range_record + decompressed_data
      #
      # Go to Complex spectral domain and window
      #
      window = np.ones(len(range_record))
      temp = 
      temp = np.fft.fft(range_record)
      #
      # Determine calibrated chirp
      # 
      calChirp = detChirpFiles(AuxDF['TX_TEMP'][_k], AuxDF['RX_TEMP'][_k])
      dechirp = np.conj(calChirp)
      
      presum_record = np.fft.
#  EDR = rangeCompression(edrname, filtering=fil)
  #
  # Make window for azimuth FFT
  #
#  window_azi = 0.5*(1.0-np.cos(2.0*np.pi*np.arange(fft_length)/(sum_length-1)))
  
#  np.save('EDR.npy', EDR)
#  print('Finished with science telemetry file...')
  return

if __name__ == '__main__':
  main()
