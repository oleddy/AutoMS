from AutoMS.automs import AutoMS
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help = 'input mzML file path', required = True)
    parser.add_argument('-o', help = 'output file path', required = True)
    parser.add_argument('-m', help = 'minimum intensity threshold', required = False, default = 5000)

    args = parser.parse_args()

    peaks = AutoMS(args.i, min_intensity = min_intensity)
    peaks.to_csv(args.o)