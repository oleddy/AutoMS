# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:15:34 2022

@author: DELL

Modified by Owen Leddy
"""

import numpy as np

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import rpy2.robjects.pandas2ri as pandas2ri

from AutoMS import peakeval

'''
To get rpy2 to work, you may need to add the following lines to your .bash_profile:
export DYLD_LIBRARY_PATH="/Library/Frameworks/R.framework/Libraries:$DYLD_LIBRARY_PATH"
export PKG_CONFIG_PATH="/usr/X11/lib/pkgconfig:$PKG_CONFIG_PATH"
'''

rcodes  =  '''
if (!'xcms' %in% installed.packages()){
  suppressMessages(install.packages('BiocManager'))
  suppressMessages(BiocManager::install('xcms'))
}
library(xcms)

getEIC <- function(file, peaks, ppm){
  tol = ppm/1000000.0
  raw_data <- readMSData(files = file, mode = "onDisk")
  rtr <- cbind(peaks$rt1, peaks$rt2)
  mzr <- cbind(peaks$mz - (peaks$mz*tol), peaks$mz + (peaks$mz*tol))
  chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
  chr_raw <- chr_raw@.Data
  output <- lapply(chr_raw, function(s) {
      if (sum(!is.na(s@intensity)) < 2){
        cbind(s@rtime, mean(s@mz), 0.0)
      }else{
        cbind(s@rtime, mean(s@mz), approx(s@intensity, n=length(s@intensity))$y)
      }
    })
  return(output)
}
'''

def AutoMS_External(file, peaks, length=14, params=(8.5101, 1.6113, 0.1950), min_width = 6, model_dir = 'model/denoising_autoencoder.pkl', ppm = 40):
    """
        1. Install R >= 3.4.1 and R <= 4.1.1
        2. Set R_HOME environment variable
        3. Install XCMS in R
    """
    numpy2ri.activate()
    pandas2ri.activate()
    
    robjects.r(rcodes)
    getEIC = robjects.globalenv['getEIC']
    
    pics_xcms = getEIC(file, peaks, ppm)
    pics_xcms = [np.array(x) for x in pics_xcms]
    
    pics_label = []
    for i, pic in enumerate(pics_xcms):
        rt, mz, intensity = peaks.loc[i, ['rt', 'mz', 'intensity']]
        label = '{}_{}_{}'.format(rt, mz, intensity)
        pics_label.append(label)
    peaks['pic_label'] = pics_label
    pics_xcms = dict(zip(pics_label, pics_xcms))
        
    scores, mspd_snrs, _, _, _ = peakeval.evaluate_peaks(peaks, pics_xcms, length=length, 
                                                                     params=params, min_width = min_width, 
                                                                     cal_snr=True, model_dir = model_dir)
    scores = np.array(scores)
    scores[scores < 0] = 0
    peaks['score'] = scores
    
    return peaks, mspd_snrs



if __name__ == '__main__':
    
    import pandas as pd
    
    file = 'data/600mix_pos.mzML'
    peaks = pd.read_csv('data/xcms_mzmine_input_demo.csv')
    peaks = AutoMS_External(file, peaks)
    
    