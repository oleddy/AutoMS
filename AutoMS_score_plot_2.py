from matplotlib import pyplot as plt
import pandas as pd
import argparse
from numpy import logical_and, logical_not

pos_ctrl_mz = ["662.312942714544", "522.7688970730040"]
pos_ctrl_indices = [1515, 2081]

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'autoMS scoring output', required = True)
    parser.add_argument('-c', help = 'mock psm peak matches', required = True)
    parser.add_argument('-m', help = 'infected psm peak matches', required = True)
   
    args = parser.parse_args()

    result = pd.read_csv(args.i)

    #import dataframes with column indicating whether a matching psm was found for mock or infected sample
    mock_matches = pd.read_csv(args.c)
    inf_matches = pd.read_csv(args.m)

    #add this information to the AutoMS sccoring results 
    result['mock_match'] = mock_matches['matching_psm']
    result['mtb_match'] = inf_matches['matching_psm']


    pos_ctrl_points = result.iloc[pos_ctrl_indices]

    mock_match_only = result.loc[logical_and(result['mock_match'], logical_not(result['mtb_match']))]
    mtb_match_only = result.loc[logical_and(result['mtb_match'], logical_not(result['mock_match']))]
    both_match = result.loc[logical_and(result['mtb_match'], result['mock_match'])]
    neither_match = result.loc[logical_and(logical_not(result['mtb_match']), logical_not(result['mock_match']))]


    print(pos_ctrl_points.shape[0])

    colors = ['gray', 'red', 'blue', 'purple']
    dfs = [neither_match, mock_match_only, mtb_match_only, both_match]

    plt.figure()
    for color, df in zip(colors, dfs):
        plt.plot(df['score_inf'], df['score_mock'], 'o', color = color)
    plt.plot(pos_ctrl_points['score_inf'], pos_ctrl_points['score_mock'], 'o', color = 'green')
    plt.legend(['neither match', 'mock match only', 'Mtb match only', 'both match', 'Known Mtb peptides'])
    plt.xlabel('Infected score')
    plt.ylabel('Mock score')
    plt.savefig('score_plot.jpeg')

    plt.figure()
    for color, df in zip(colors, dfs):
        plt.plot(df['snr_inf'], df['snr_mock'], 'o', color = color)
    plt.plot(pos_ctrl_points['snr_inf'], pos_ctrl_points['snr_mock'], 'o', color = 'green')
    plt.legend(['neither match', 'mock match only', 'Mtb match only', 'both match', 'Known Mtb peptides'])
    plt.xlabel('Infected SNR')
    plt.ylabel('Mock SNR')
    plt.savefig('snr_plot.jpeg')
