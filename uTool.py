from pyopenms import *

#Get data from input file
exp = MSExperiment()
MzMLFile().load("tiny.pwiz.1.1.mzML", exp)

#Convert the string representation of the peptide into an amino acid sequence object
#Hardcoded for now, get input later?
ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
peptide = AASequence.fromString(ubiquitin)

tsg = TheoreticalSpectrumGenerator()
spec = MSSpectrum()
p = Param()
p.setValue(b"add_metainfo", b"true", "")
tsg.setParameters(p)
tsg.getSpectrum(spec, peptide, 1, 2)
print("Spectrum has", spec.size(), "peaks.")
print(spec.get_peaks())
#for ion, peak in zip(spec.getStringDataArrays()[0], spec):
#   print(ion, peak.getMZ())