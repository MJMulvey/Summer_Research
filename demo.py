#TODO Why do I have to import both? (because of pyopenms.msexperiment() and msexperiment?)
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
def get_difference(bound, unbound):
    bound_peaks = {}
    unbound_peaks = {}
    difference_peaks = {}
    #TODO Resolution - does get_peaks provide enough? Too much?
    #TODO How should I account for small variations / noise? e.g. Two peaks in result and original with the same intensity but a tiny difference in mz (which may be the same molecule in reality)
    for mz, intensity in bound.items():
        if mz in bound_peaks:
            bound_peaks[mz] += intensity
        else:
            bound_peaks[mz] = intensity
    for mz, intensity in unbound.items():
        if mz in unbound_peaks:
            unbound_peaks[mz] += intensity
        else:
            unbound_peaks[mz] = intensity
    bound_max = max(bound_peaks.values())
    unbound_max = max(unbound_peaks.values())
    for mz, i in bound_peaks.items():
        difference_peaks[mz] = i / bound_max
    for mz, i in unbound_peaks.items():
        if mz in difference_peaks:
            difference_peaks[mz] -= i / unbound_max
        else:
            difference_peaks[mz] = - i / unbound_max
    return(difference_peaks)

#Remove line breaks from plain test data
def remove_breaks(lines):
    for index in range(len(lines)):
        line = lines[index]
        if line[-1] == "\n":
            lines[index] = line[0:-1]
    return(lines)

#Open the unbound file and read the lines
#unbound_file = open("Ub_ISO779_1.xy", "r")
unbound_file = open("ub_5_1.xy", "r")
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
#bound_file = open("Ub_C_ISO779_1.xy", "r")
bound_file = open("t_2.xy", "r")
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

#Convert the string representation of Ubiquitin into an amino acid sequence object
ubiquitin = AASequence.fromString("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
tsg = TheoreticalSpectrumGenerator()
spectrum = MSSpectrum()
#Initialise parameters object and set the values for parameters
parameters = Param()
parameters.setValue(b"add_isotopes", b"true", "")
parameters.setValue(b"add_losses", b"true", "")
parameters.setValue(b"add_b_ions", b"true", "")
parameters.setValue(b"add_y_ions", b"true", "")
parameters.setValue(b"add_a_ions", b"true", "")
parameters.setValue(b"add_c_ions", b"true", "")
parameters.setValue(b"add_x_ions", b"true", "")
parameters.setValue(b"add_z_ions", b"true", "")
parameters.setValue(b"add_precursor_peaks", b"true", "")
parameters.setValue(b"add_all_precursor_charges", b"true", "")
parameters.setValue(b"add_first_prefix_ion", b"true", "")
parameters.setValue(b"add_metainfo", b"true", "")
tsg.setParameters(parameters)
#Generate the theoretical spectrum of Ubiquitin
tsg.getSpectrum(spectrum, ubiquitin, 1, 2)

from scipy.spatial.distance import euclidean
from fastdtw import fastdtw

difference = []
for mz, i in binding_effect.items():
    if i < 0:
        difference.append([mz,-i])

theoretical = []
for i in range(len(spectrum.get_peaks()[0])):
    current_mz = spectrum.get_peaks()[0][i]
    current_intensity = spectrum.get_peaks()[1][i]
    theoretical.append([current_mz, current_intensity])

from operator import itemgetter

length = len(theoretical)
filtered_difference = sorted(difference, key=itemgetter(1), reverse=True)
filtered_difference = filtered_difference[:length]

#import matplotlib.pyplot as plt
#for i in range(0,len(filtered_difference),1):
#    item = filtered_difference[i]
#    plt.vlines(item[0],0,item[1])
#    print(i)
#plt.show()

#print("Diff", len(difference))
#print("Theo", len(theoretical))
#print("Filt", len(filtered_difference))

distance, path = fastdtw(difference, theoretical)
norm_distance = distance / length
print("Dist", distance)
print("Norm dist", norm_distance)

#input_map = MSExperiment()
#MzMLFile().load("bound.mzML", input_map)
#input_map.updateRanges()