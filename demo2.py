from pyopenms import *

exp = MSExperiment()
MzMLFile().load("velos005614.mzML", exp)
spectra = exp.getSpectrum(1)
print(spectra.get_peaks())