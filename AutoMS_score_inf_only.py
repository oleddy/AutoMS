import pandas as pd
import argparse

from os.path import join
import os
import sys

script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
script_directory = script_directory.replace('/scripts', '') #for some reason the above gives script_directory as AutoMS/scripts instead of AutoMS even though that's where the script file is

from AutoMS.automs_external import AutoMS_External

from math import ceil
from functools import partial
import multiprocessing as mp

mp.set_start_method('fork') #prevents a file not found error (idk why)

#breaks down a large dataframe of MS1 peaks into smaller chunks to avoid running out of memory by running AutoMS on all at once
def chunk(df, chunksize = 4096):
    return [df.iloc[i*chunksize:(i+1)*chunksize].copy().reset_index(drop = True) for i in range(ceil(len(df)/chunksize))] #need to reset index because AutoMS assumes first index of dataframe is 0

#function to score each chunk of the larger dataframe
def score_chunk(peaks, inf_mzml_file, ppm = 40, length = 40): 
    model_dir = join(script_directory, 'model/denoising_autoencoder.pkl')
    inf_result, inf_snrs = AutoMS_External(inf_mzml_file, peaks, min_width = 0, length = float(length), model_dir = model_dir, ppm = int(ppm))
    result = inf_result.rename({'score' : 'score_inf'}, axis = 1)
    result['snr_inf'] = inf_snrs #add signal/noise ratio as an additional column in results dataframe
    return result
    
#score full dataframe by breaking it into smaller chunks and scoreing the chunks, then rejoining them
def AutoMS_score(unmatched_peaks, inf_mzml_file, outfile, ppm = 40, length = 40, n_cores = mp.cpu_count() - 2, chunksize = 4096):
    peaks = pd.read_csv(unmatched_peaks)
    chunks = chunk(peaks, chunksize)
    score_partial = partial(score_chunk, inf_mzml_file = inf_mzml_file, ppm = ppm, length = length)   
    with mp.Pool(n_cores) as pool:
        results = pool.map(score_partial, chunks)
    # results = [score_partial(chunk) for chunk in chunks] #version without multiprocessing
    result = pd.concat(results, ignore_index = True)
    result.to_csv(outfile, index = False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', help = 'input (unpaired peaks file)', required = True)
    parser.add_argument('-i', help = 'infected MZML file', required = True)
    parser.add_argument('-o', help = 'output path', required = True)
    parser.add_argument('-w', help = 'mass window for XIC, in ppm', required = False, default = 40, type = float)
    parser.add_argument('-l', help = 'XIC time window', required = False, default = 40, type = float)
    parser.add_argument('-c', help = 'chunk size (adjust to avoid hitting memory limits)', type = int, required = False, default = 4096)
    parser.add_argument('-n', help = 'number of cores', required = False, type = int, default = mp.cpu_count() - 2)

    args = parser.parse_args()
    
    results = AutoMS_score(args.u, args.i, args.o, ppm = args.w, length = args.l, chunksize = args.c)
