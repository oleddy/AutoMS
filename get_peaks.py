import argparse
import pandas as pd
import numpy as np
import multiprocessing as mp
from math import ceil
from functools import partial

'''
This script takes the peaks output from AutoMS (a peak 
calling denoising autoencoder), and returns a tab-delimited 
list of peaks in suitable format for input to chromatographic
alignment algorithm, as well as the filtered peak list from 
AutoMS with charge states added
'''

#average mass gain with addition of one neutron, in Da
avg_neutron = 1.002

#for a given base peak, test whether a given other peak is an isotopic peak of the base peak, for a given charge state
def is_isotopic(base_peak, test_peak, z, mz_tolerance = 0.01, RT_tolerance = 60.): #NOTE: how tight should these tolerances be?
    #check if the offset between the base peak and some test peak is equal to 1/z, within some tolerance
    mz_match = abs((test_peak['mz'] - base_peak['mz']) - (avg_neutron/z)) <= mz_tolerance 

    #check if the estimated RT of the two peaks is within some tolerance
    RT_match = abs(base_peak['rt'] - test_peak['rt']) <= RT_tolerance 

    #check if the intensity of the test peak is lower than that of the base peak (NOTE: Should this be required?? May be some cases where 13C peak is higher than base)
    lower_intensity = base_peak['intensity'] > test_peak['intensity'] 
    
    if mz_match and RT_match and lower_intensity:
        return True

#size of window (in m/z units) to examine for isotopic peaks
isotope_window = 3. 

#function that takes a list of peaks output by AutoMS and groups them into sets of isotopic peaks to assign charge states
def group_isotopic_peaks(range, charges):

    #unpack the argument range (specifying index range) into beginning and end index
    begin, end = range

    #retrieve the appropriate set of peaks from the peaks dataframe, set as global variable by multiprocess_df() (see below)
    peaks = df_mp.iloc[begin:end].reset_index(drop = True)

    #make sure we test the various charge states from highest to lowest, so that if we have a +4 charge state, for example, we don't assign it as +2
    charges.sort()
    charges.reverse()

    #initialize a dataframe to hold the base peaks with charge states assigned
    grouped_peaks = pd.DataFrame(columns = list(peaks.columns) + ['charge']) 

    #transpose so we can use the pop() method to pop off the next peak and remove it from the dataframe
    peaks_T = peaks.T

    #we're going to pop off each peak as we define its isotopic series (or fail to find any isotopic peaks)
    #so we're done when no peaks are left
    while not peaks_T.empty:

        #pop off next lowest m/z peak
        current_peak = peaks_T.pop(0)

        #iterate over charges (from greatest to least) until we find an isotopic series or hit z = 1
        for z in charges: 
            isotopic_set = [current_peak]

            #examine only peaks within a window extending from the mz of the current peak up to a finite window (width defined by isotope_window)
            mz_close = np.logical_and(peaks_T.loc['mz'] > current_peak['mz'], peaks_T.loc['mz'] < current_peak['mz'] + isotope_window)
            close_peaks = peaks_T.loc[:,mz_close]
            
            #iterate over peaks within the appropriate range to find next isotopic peak. If we hit the end of the window without finding one, we're done
            for i, peak2 in close_peaks.iteritems():
                if is_isotopic(current_peak, peak2, z):
                    next_peak = peaks_T.pop(i)
                    isotopic_set.append(next_peak)
                    current_peak = next_peak

            #if we found at least one isotopic peak, define this is a group with a known charge state and add it to the list of grouped peaks
            if len(isotopic_set) > 1:
                grouped_peak = isotopic_set[0] #copy base peak
                grouped_peak['charge'] = z #add charge information
                grouped_peaks = pd.concat([grouped_peaks, pd.DataFrame(grouped_peak).T], axis = 0) #append to dataframe of peak groups
                break #if we've found a series of isotopic peaks, stop iterating over charge states for that particular base peak

        charged_peak = current_peak
        charged_peak['charge'] = 1 #assign charge state 1 by default if no isotopic peaks found (NOTE: is this the best solution?? is it better to assign 0 maybe, or will SIMA accept NA here?)
        grouped_peaks = pd.concat([grouped_peaks, pd.DataFrame(charged_peak).T], axis = 0)

        if not peaks_T.empty:
            #transpose, reset indices, and transpose back
            peaks_T = peaks_T.T.reset_index(drop = True).T

    return grouped_peaks

#after adding charge states, converts the output dataframe from AutoMS to a dataframe suitable for SIMA format (when output as a tab-delimited file)
#path to the original mzML file is needed to retrieve charge of each precursor
def sima_format(peaks):
    peaks_sima = peaks[['charge', 'intensity', 'rt']].copy()
    peaks_sima['mass'] = peaks['mz']*peaks['charge']
    peaks_sima = peaks_sima[['mass', 'charge', 'intensity', 'rt']] #reorder columns
    return peaks_sima

#helper function for multiprocess_df() below. initialize the dataframe as a global variable accessible by multiprocessing child processes
#so we don't have to turn it into a bytestream and put it through a pipe to each child process, which takes a long time if the df is large
def multiprocess_df_initialize(_df):
    global df_mp
    df_mp = _df

#split a dataframe into chunks, process each chunk (using a function assumed to return a pandas DataFrame) and reassemble
def multiprocess_df(df, fxn, **kwargs):
    #use maximum allowed number of cores (reserving one for parent process and one for non-python tasks)
    n_cores = mp.cpu_count() - 2

    #number of rows per chunk (rounded up)
    n_per_chunk = ceil(df.shape[0]/n_cores)

    #split df into a list of chunks
    chunks = [(i*n_per_chunk,(i+1)*n_per_chunk) for i in range(n_cores)]

    #mapping with multiprocessing requires single-function argument, so pre-apply other arguments besides chunk itself
    fxn_partial = partial(fxn, **kwargs)

    # multiprocess chunks using function with pre-applied arguments and with the dataframe initialized as a global variable via multiprocess_df_initialize()
    with mp.Pool(n_cores, initializer = multiprocess_df_initialize, initargs = (df,)) as pool:
        chunks_processed = pool.map(fxn_partial, chunks)
    
    #reassemble chunks and return (ignore index since fxn may have reset the indices of each chunk to begin at 0)
    return pd.concat(chunks_processed, ignore_index = True)

charges = [4, 3, 2] #charge states for isotopic peak grouping. NOTE: do we only care about charge states that MHC peptides are likely to have? 
#we can omit z=1 to save time since we assign it by default if we find no isotopic peaks with the other z values

#put it all together -- read in AutoMS output results and output the filtered peaks with charge states assigned in both the original format and SIMA format
def get_peaks(file, score_threshold, output_path_raw, output_path_sima):
    #load autoMS output
    peaks = pd.read_csv(file)
    
    #sort to make sure we go from lowest m/z to highest (avoids picking isotopic peaks before base peak)
    peaks = peaks.sort_values(['mz']) 
    peaks = peaks.reset_index(drop = True)

    #take only peaks over score threshold
    peaks_thresh = peaks.loc[peaks['score'] >= score_threshold]
    peaks_thresh = peaks_thresh.reset_index(drop = True)

    #add charge states to peaks dataframe
    peaks_charge = multiprocess_df(peaks_thresh, group_isotopic_peaks, charges = charges)

    #save peaks dataframe
    peaks_charge.to_csv(output_path_raw) 

    #convert to sima format
    peaks_sima = sima_format(peaks_charge)

    #save sima formatted peaks
    peaks_sima.to_csv(output_path_sima, sep = '\t', index = False)
    

#if calling from command line, read in arguments
if __name__ == '__main__':

    #parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'input autoMS output file', required = True)
    parser.add_argument('-t', help = 'threshold quality score', required = False, default = 1.)
    parser.add_argument('-o', help = 'raw output path', required = True)
    parser.add_argument('-s', help = 'SIMA output path', required = True)

    args = parser.parse_args()

    get_peaks(args.i, float(args.t), args.o, args.s)
    