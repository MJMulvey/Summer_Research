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
input_file = open("\xy\ub\ub_5_2.xy", "r")
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

import matplotlib.pyplot as plt
x = list(pairs.keys())
y = list(pairs.values())
for i in range(0,len(x),500):
    plt.vlines(x[i],0,y[i])
    print(i)
plt.show()