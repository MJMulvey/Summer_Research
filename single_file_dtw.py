from pyopenms import *
import pyopenms

#Returns a dictionary of peaks obtained from taking the peaks of the original spectrum away from the peaks of the result spectrum
def normalise(mzs, intensities):
    pairs = {}
    for i in range(len(mzs)):
        current_mz = mzs[i]
        current_intensity = intensities[i]
        if current_mz in pairs:
            pairs[current_mz] += current_intensity
        else:
            pairs[current_mz] = current_intensity
    maximum = max(pairs.values())
    for current_mz, current_intensity in pairs.items():
        pairs[current_mz] = current_intensity / maximum
    return(pairs)

#Remove line breaks from plain test data
def remove_breaks(lines):
    for index in range(len(lines)):
        line = lines[index]
        if line[-1] == "\n":
            lines[index] = line[0:-1]
    return(lines)

#Open the file and read the lines
input_file = open("demo_4.xy", "r")
input_lines = input_file.readlines()
input_file.close()
input_lines = remove_breaks(input_lines)

mz = []
intensity =[]
#Process each line of the text file and store it in arrays
for line in input_lines:
    current_mz, current_intensity = line.split()
    mz.append(float(current_mz))
    intensity.append(float(current_intensity))

pairs = normalise(mz, intensity)

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
for mz, i in pairs.items():
    difference.append([mz,i])

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

#print("Theo", len(theoretical))
#print("Filt", len(filtered_difference))

len_fd = len(filtered_difference)

distance, path = fastdtw(filtered_difference, theoretical)
norm_distance = distance / (len_t + len_fd)
print("Dist", distance)
print("Norm dist", norm_distance)
print("Board =", path[-1])

#import matplotlib.pyplot as plt
#for i in range(0,len(theoretical),1):
#    item = theoretical[i]
#    plt.vlines(item[0],0,item[1])
#    print(i)
#plt.show()