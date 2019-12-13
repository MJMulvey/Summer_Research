#TODO Why do I have to import both?
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

#Remove line breaks from plain test data
def remove_breaks(lines):
    for index in range(len(lines)):
        line = lines[index]
        if line[-1] == "\n":
            lines[index] = line[0:-1]
    return(lines)

#Open the unbound file and read the lines
unbound_file = open("Ub_ISO779_1.xy", "r")
unbound_lines = unbound_file.readlines()
unbound_file.close()
unbound_lines = remove_breaks(unbound_lines)
#Create an experiment to store the unbound data
unbound_exp = pyopenms.MSExperiment()
unbound_spectrum = MSSpectrum()
unbound_mz = []
unbound_intensity =[]
#Process each line of the text file and store it in arrays
for line in unbound_lines:
    current_mz, current_intensity = line.split()
    unbound_mz.append(float(current_mz))
    unbound_intensity.append(float(current_intensity))
#Update the experiment with the data, then store it in a file
unbound_spectrum.set_peaks([unbound_mz, unbound_intensity])
unbound_exp.setSpectra([unbound_spectrum])
pyopenms.MzMLFile().store("unbound.mzML", unbound_exp)

#Open the bound file and read the lines
bound_file = open("Ub_C_ISO779_1.xy", "r")
bound_lines = bound_file.readlines()
bound_file.close()
bound_lines = remove_breaks(bound_lines)
#Create an experiment to store the bound data
bound_exp = pyopenms.MSExperiment()
bound_spectrum = MSSpectrum()
bound_mz = []
bound_intensity =[]
#Process each line of the text file and store it in arrays
for line in bound_lines:
    current_mz, current_intensity = line.split()
    bound_mz.append(float(current_mz))
    bound_intensity.append(float(current_intensity))
#Update the experiment with the data, then store it in a file
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

#Calculate the effect of binding on the unbound spectrum
binding_effect = get_difference(bound_spectrum, unbound_spectrum)

output_file = open("output.txt", "w")
output_string = ""
for mz, i in binding_effect.items():
    output_string += str(mz) + "\t" + str(i) + "\n"
output_file.write(output_string)
output_file.close()