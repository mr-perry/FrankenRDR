#
# Import necessary libraries
#
import os
import struct
import bitstring
import numpy as np


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
