import pandas as pd
import argparse
from AutoMS.automs_external import AutoMS_External

min_width = 0
length = 40
ppm = 10

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', help = 'input (unpaired peaks file)', required = True)
    parser.add_argument('-i', help = 'infected MZML file', required = True)
    parser.add_argument('-c', help = 'control (mock infected) MZML file', required = True)
    parser.add_argument('-o', help = 'output path', required = True)

    args = parser.parse_args()

    peaks = pd.read_csv(args.u)
    # model_dir = '../AutoMS/model/denoising_autoencoder.pkl'
    model_dir = 'model/denoising_autoencoder.pkl'
    inf_result, inf_snrs = AutoMS_External(args.i, peaks, min_width = min_width, length = length, model_dir = model_dir, ppm = ppm)
    result = inf_result.rename({'score' : 'score_inf'}, axis = 1)
    result['snr_inf'] = inf_snrs #add signal/noise ratio as an additional column in results dataframe

    mock_result, mock_snrs = AutoMS_External(args.c, peaks, min_width = min_width, length = length, model_dir = model_dir, ppm = ppm)
    result['snr_mock'] = mock_snrs #add signal/noise ratio as an additional column in results dataframe
    result['score_mock'] = mock_result['score']

    result.to_csv(args.o, index = False)
