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
            difference_peaks[mz] = -i / unbound_max
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
unbound_file = open("ub_1_1.xy", "r")
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
    #Accuracy of mz values is standardised at 1dp as this is the lowest accuracy encountered
    unbound_mz.append(round(float(current_mz),1))
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
    #Accuracy of mz values is standardised at 1dp as this is the lowest accuracy encountered
    bound_mz.append(round(float(current_mz),1))
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
#parameters.setValue(b"add_losses", b"true", "")
parameters.setValue(b"add_b_ions", b"true", "")
parameters.setValue(b"add_y_ions", b"true", "")
parameters.setValue(b"add_a_ions", b"true", "")
parameters.setValue(b"add_c_ions", b"true", "")
parameters.setValue(b"add_x_ions", b"true", "")
parameters.setValue(b"add_z_ions", b"true", "")
#parameters.setValue(b"add_precursor_peaks", b"true", "")
#parameters.setValue(b"add_all_precursor_charges", b"true", "")
#parameters.setValue(b"add_first_prefix_ion", b"true", "")
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

len_t = len(theoretical)
filtered_difference = sorted(difference, key=itemgetter(1), reverse=True)
if len_t < len(filtered_difference):
    filtered_difference = filtered_difference[:len_t]
#filtered_difference = filtered_difference[:1000]

len_fd = len(filtered_difference)

distance, path = fastdtw(filtered_difference, theoretical)
#Using M + N (M = the length of one side, and N = the length of the other) because it's proportional to the length of the diagonal across the grid
norm_distance = distance / (len_t + len_fd)
#print("Dist", distance)
#print("Norm dist", norm_distance)

#Tried different numbers of significant peaks (s). In one sample the following s / (lowest % intensity captured) pairs were obtained
#(10,65%),(25,29%),(50,12%),(100,2.9%),(200,1.0%),(300,0.90%),(500,0.78%),(1000,0.51%)
significant_peaks = filtered_difference[:100]

#Values in matches are indicies
#matches = {}
#for i in range(len(significant_peaks)):
#    matches[i] = []
#for link in path:
#    if link[0] in matches:
#        matches[link[0]].append(link[1])

#rounded_theoretical = {}
#for i in range(len(theoretical)):
#    peak = theoretical[i]
#    rounded = round(peak[0],1)
#    if rounded not in rounded_theoretical.keys():
#        rounded_theoretical[rounded] = peak[1]
    #Duplicate theoretical peaks with the same rounded mz are discarded (as we're unable to tell which of them the experimental peak matches to)
    #The peak with the largest intensity is retained
#    else:
#        rounded_theoretical[rounded] = max(peak[1],rounded_theoretical[rounded])

def find_closest(peak):
    if peak[0] < theoretical[0][0]:
        return(theoretical[0])
    elif peak[0] > theoretical[-1][0]:
        return(theoretical[-1])
    else:
        left = 0
        right = 0
        index = 0
        while left == 0 :
            current = theoretical[index]
            if current[0] == peak[0]:
                return(current)
            elif current[0] > peak[0]:
                right = current
                left = theoretical[index - 1]
            index += 1
        left_closer = peak[0] - left[0] < right[0] - peak[0]
        if left_closer:
            return(left)
        else:
            return(right)

matched_peaks = []
for peak in significant_peaks:
    theo_peak = find_closest(peak)
    matched_peaks.append([peak, theo_peak])


matching_significance = {}
for match in matched_peaks:
    charge = match[1][0]
    intensity = match[0][1]
    if charge not in matching_significance.keys():
        matching_significance[charge] = intensity
    else:
        matching_significance[charge] += intensity

fragments = []

for ion, peak in zip(spectrum.getStringDataArrays()[0], spectrum):
    for peak_mz, peak_sig in matching_significance.items():
        if peak.getMZ() == peak_mz:
            fragments.append([ion, peak_mz, peak_sig])

fragments = sorted(fragments, key=itemgetter(2), reverse=True)

print("Ion\t\tm/z\t\t\tRelative Significance")
for fragment in fragments:
    if len(fragment[0]) < 5:
        print(fragment[0], "\t\t", fragment[1], "\t\t", fragment[2], sep="")
    else:
        print(fragment[0], "\t", fragment[1], "\t\t", fragment[2], sep="")