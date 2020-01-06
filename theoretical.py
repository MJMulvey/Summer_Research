from pyopenms import *

#Convert the string representation of Ubiquitin into an amino acid sequence object
ubiquitin = AASequence.fromString("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
tsg = TheoreticalSpectrumGenerator()
spectrum = MSSpectrum()
#Initialise parameters object and set the values for parameters
parameters = Param()
parameters.setValue(b"add_b_ions", b"true", "")
parameters.setValue(b"add_y_ions", b"true", "")
parameters.setValue(b"add_a_ions", b"true", "")
parameters.setValue(b"add_c_ions", b"true", "")
parameters.setValue(b"add_x_ions", b"true", "")
parameters.setValue(b"add_z_ions", b"true", "")
parameters.setValue(b"add_precursor_peaks", b"true", "")
parameters.setValue(b"add_metainfo", b"true", "")
tsg.setParameters(parameters)
tsg.getSpectrum(spectrum, ubiquitin, 1, 2)
print("Spectrum has", spectrum.size(), "peaks.")
#print(spectrum.get_peaks())
#for ion, peak in zip(spectrum.getStringDataArrays()[0], spectrum):
#   print(ion, peak.getMZ())