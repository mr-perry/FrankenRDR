#
# Import necessary libraries
#
import numpy as np
import bitstring


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
  decoded_data = np.zeros(3600, int)
  #
  # Step 4: Decode the science data based on the bit resolution
  #
  # Break Byte stream into Bit Stream
  # This isn't necessary for the 8-bit resolution, but to keep the code
  # clean, split the bytes up regardless of resolution
  #
  b = bitstring.BitArray(echoes)
  for _j in range(0, len(b), bps):
    decoded_data[int(_j/bps)] = b[_j:_j+bps].int
  return decoded_data, ancil


def decompressSciData(data, compression, presum, bps):
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
    L = np.ceil(np.log2(int(presum)))
    R = bps
    S = L - R + 8
    N = presum
    decompressed_data = data * (np.power(2, S)/N)
    return decompressed_data
  elif compression == True:#dynamic scaling
    print('This is not yet available.')
    sys.exit()
