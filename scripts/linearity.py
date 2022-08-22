# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 10:01:39 2022

@author: DELL
"""

import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm

from AutoMS import hpic
from AutoMS import peakeval

path = 'E:/Data/Quantification Dataset for KPIC2"'
files = os.listdir(path)

output = {}
for f in tqdm(files):
    file = path + '/' + f
    peaks, pics = hpic.hpic(file, min_intensity=500, min_snr=1, mass_inv=1, rt_inv=30)
    scores, _, _, _, _, _ = peakeval.evaluate_peaks(peaks, pics, cal_snr=False, use_cnn=False)
    peaks['scores'] = scores
    output[f] = peaks