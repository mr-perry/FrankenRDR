#
# Import necessary libraries
#
import numpy as np

def makeWindow(win_type, len=2048):
  """
    Function for writing a windowing function
    Inputs:
      win_type: 0 -- Uniform weighting
                1 -- Cosine Bell (not working)
                2 -- Hann Window
                3 -- Hamming Window
  *** NOTE *** This will eventually include a generalized Hamming function 
  """
  #
  # Construct the window function
  #
  win_str ='Window function applied to data: '
  if win_type == 0:
    window = np.ones(2048)                                      # Uniform weighting
    win_str += 'Uniform weighting'
  elif win_type == 1:
    #
    # Cosine bell
    window = np.ones(2048)                                      # Uniform weighting
    win_str += 'Uniform weighting'
  elif win_type == 3:
    window = 0.5 * (1.0 - np.cos(2.0*np.pi*(np.arange(0,2048,1) - 2048)/2047))
    win_str += 'Hann Window'
  else:
    writeLog(_log, 'WARNING: No windowing function selected!')
  return window, win_str
