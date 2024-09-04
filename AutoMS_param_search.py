import pandas as pd
import argparse
from AutoMS.automs_external import AutoMS_External

min_width = 0 #setting any value >0 for minimum peak width appears to exclude some peaks we would want to capture, since some peaks are only present in 1-2 MS1 scans
lengths = [10, 20, 40, 60, 100]
ppms = [10, 20, 40]

'''
Run AutoMS on a dataset with a range of parameter combinations (mass tolerance and XIC length) and save the results. 
'''

#TODO: add multiprocessing
def add_scores(length, ppm, mzml_file, peaks, peaks_output):
    model_dir = 'model/denoising_autoencoder.pkl'
    result, snrs = AutoMS_External(mzml_file, peaks, min_width = min_width, length = length, model_dir = model_dir, ppm = ppm)
    peaks_output['score_len' + str(length) + '_' + str(ppm) + 'ppm'] = result['score']
    peaks_output['snr_len' + str(length) + '_' + str(ppm) + 'ppm'] = snrs

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', help = 'input (unpaired peaks file)', required = True)
    parser.add_argument('-i', help = 'infected MZML file', required = True)
    # parser.add_argument('-c', help = 'control (mock infected) MZML file', required = True)
    parser.add_argument('-o', help = 'output path', required = True)

    args = parser.parse_args()

    peaks = pd.read_csv(args.u)
    peaks_output = peaks.copy() #make a copy of the dataframe so we can add the scores here without cluttering the dataframe that gets passed to AutoMS_external

    for ppm in ppms:
        for length in lengths:
            add_scores(length, ppm, args.i, peaks, peaks_output)
            print(peaks_output.head())
    # add_scores_partial = partial(add_scores, mzml_file = args.i, peaks = peaks, peaks_output = peaks_output)
    
    # with mp.Pool(4) as pool:
    #     _ = pool.map(add_scores_partial, param_pairs)

    peaks_output.to_csv(args.o, index = False)
