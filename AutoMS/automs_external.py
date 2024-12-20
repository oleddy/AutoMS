# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 17:15:34 2022

@author: DELL

Modified by Owen Leddy
"""

import numpy as np
import pandas as pd

from AutoMS import peakeval

import pyopenms as oms

def getXIC(file, peaks, ppm, length):
    run = oms.MSExperiment()
    oms.MzMLFile().load(file, run)
    peak_outputs = [pd.DataFrame(columns = ['rt', 'mz', 'intensity']) for i in range(len(peaks))]

    for spectrum in run:
        rt = spectrum.getRT()
        for i, peak in peaks.iterrows():
            if (rt <= peak['rt'] + length) and (rt >= peak['rt'] - length):
                tolerance = peak['mz']*(ppm/1e6)
                index = spectrum.findHighestInWindow(peak['mz'], tolerance, tolerance)
                if index == -1:
                    intensity = 0.
                else:
                    intensity = spectrum[index].getIntensity()
                new_row = pd.DataFrame({'rt' : [rt], 'mz' : [peak['mz']], 'intensity' : [intensity]})
                peak_outputs[i] = pd.concat([peak_outputs[i], new_row])
    return peak_outputs    

def AutoMS_External(file, peaks, length=14, params=(8.5101, 1.6113, 0.1950), min_width = 6, model_dir = 'model/denoising_autoencoder.pkl', ppm = 40):
    pics_xcms = getXIC(file, peaks, ppm, length)
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
    
    