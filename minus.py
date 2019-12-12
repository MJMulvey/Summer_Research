from pyopenms import *
import pyopenms

#Combines several spectra (separated by time) into a single spectrum for that experiment
def sum_spectra(spectra):
    summed_peaks = {}
    for spectrum in spectra:
        mz, intensity = spectrum.get_peaks()
        for i in range(len(mz)):
            current_mz = mz[i]
            current_intensity = intensity[i]
            if current_mz in summed_peaks:
                summed_peaks[current_mz] += current_intensity
            else:
                summed_peaks[current_mz] = current_intensity
    return(summed_peaks)

#Returns a dictionary of peaks obtained from taking the peaks of the original spectrum away from the peaks of the result spectrum
def get_difference(result, original):
    difference_peaks = {}
    #TODO Resolution - does get_peaks provide enough? Too much?
    #TODO How should I account for small variations / noise? e.g. Two peaks in result and original with the same intensity but a tiny difference in mz (which may be the same molecule in reality)
    for mz, intensity in result.items():
        if mz in difference_peaks:
            difference_peaks[mz] += intensity
        else:
            difference_peaks[mz] = intensity
    for mz, intensity in original.items():
        if mz in difference_peaks:
            difference_peaks[mz] -= intensity
        else:
            difference_peaks[mz] = -intensity
    return(difference_peaks)

def remove_breaks(lines):
    for index in range(len(lines)):
        line = lines[index]
        if line[-1] == "\n":
            lines[index] = line[0:-1]
    return(lines)

unbound_file = open("Ub_trial.xy", "r")
unbound_lines = unbound_file.readlines()
unbound_file.close()
unbound_lines = remove_breaks(unbound_lines)

unbound_exp = pyopenms.MSExperiment()
unbound_spectrum = MSSpectrum()
unbound_mz = []
unbound_intensity =[]

for line in unbound_lines:
    current_mz, current_intensity = line.split()
    unbound_mz.append(float(current_mz))
    unbound_intensity.append(float(current_intensity))

unbound_spectrum.set_peaks([unbound_mz, unbound_intensity])
unbound_exp.setSpectra([unbound_spectrum])
pyopenms.MzMLFile().store("unbound.mzML", unbound_exp)

bound_file = open("Ub_c_trial.xy", "r")
bound_lines = bound_file.readlines()
bound_file.close()
bound_lines = remove_breaks(bound_lines)

bound_exp = pyopenms.MSExperiment()
bound_spectrum = MSSpectrum()
bound_mz = []
bound_intensity =[]

for line in bound_lines:
    current_mz, current_intensity = line.split()
    bound_mz.append(float(current_mz))
    bound_intensity.append(float(current_intensity))

bound_spectrum.set_peaks([bound_mz, bound_intensity])
bound_exp.setSpectra([bound_spectrum])
pyopenms.MzMLFile().store("bound.mzML", bound_exp)

#Import the files for bound and unbound Ubiquitin into experiments, and get an i / mz spectrum for each experiment
bound = MSExperiment()
MzMLFile().load("bound.mzML", bound)
bound_spectrum = sum_spectra(bound.getSpectra())
unbound = MSExperiment()
MzMLFile().load("unbound.mzML", unbound)
unbound_spectrum = sum_spectra(unbound.getSpectra())

binding_effect = get_difference(bound_spectrum, unbound_spectrum)

print(binding_effect[0])