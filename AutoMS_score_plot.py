from matplotlib import pyplot as plt
import pandas as pd
import argparse

pos_ctrl_mz = ["662.312942714544", "522.7688970730040"]
pos_ctrl_indices = [1515, 2081]

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'autoMS scoring output', required = True)
   
    args = parser.parse_args()

    result = pd.read_csv(args.i)

    pos_ctrl_points = result.iloc[pos_ctrl_indices]
    print(pos_ctrl_points.shape[0])

    plt.figure()
    plt.plot(result['score_inf'], result['score_mock'], 'o')
    plt.plot(pos_ctrl_points['score_inf'], pos_ctrl_points['score_mock'], 'o')
    plt.legend(['All peptides', 'Known Mtb peptides'])
    plt.xlabel('Infected score')
    plt.ylabel('Mock score')
    plt.savefig('score_plot.jpeg')

    plt.figure()
    plt.plot(result['snr_inf'], result['snr_mock'], 'o')
    plt.plot(pos_ctrl_points['snr_inf'], pos_ctrl_points['snr_mock'], 'o')
    plt.legend(['All peptides', 'Known Mtb peptides'])
    plt.xlabel('Infected SNR')
    plt.ylabel('Mock SNR')
    plt.savefig('snr_plot.jpeg')
