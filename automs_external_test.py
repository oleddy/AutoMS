import pandas as pd
from AutoMS.automs_external import AutoMS_External

mzml_file = '../mtb_donorA_test/data/inf.mzML'
peaks_file = '../mtb_donorA_test/unpaired_peaks_automs.csv'

peaks = pd.read_csv(peaks_file)

result, snrs = AutoMS_External(mzml_file, peaks)
result.to_csv('automs_external_test.csv')
pd.DataFrame(snrs).to_csv('snr.csv')