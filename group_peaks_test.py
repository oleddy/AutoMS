from get_peaks import group_isotopic_peaks
import pandas as pd

if __name__ == '__main__':
    peaks = pd.read_csv('test peaks.csv')
    charges = [1, 2, 3, 4]
    grouped_peaks = group_isotopic_peaks(peaks, charges)
    print(grouped_peaks.head())
    