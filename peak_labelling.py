from pyopenms import *
from operator import itemgetter

#Remove line breaks from plain text data
def remove_breaks(lines):
    for index in range(len(lines)):
        line = lines[index]
        if line[-1] == "\n":
            lines[index] = line[0:-1]
    return(lines)

#Filtering should be improved, maybe based on a distance in i (related to the height of the peak?) rather that 5 L/R?
def get_maxima(peaks):
    maxima = []
    appendable = True
    for i in range(len(peaks)):
        peak = peaks[i]
        #This is a significance parameter for peaks (change it?)
        if peak[1] > 0.05 and appendable:
            maximum = True
            #Left and right bound values of being 5 away are also arbitrary parameters
            left_bound = max(0, i - 20)
            right_bound = min(len(peaks), i + 20)
            for j in range(left_bound, right_bound):
                comparison_peak = peaks[j]
                if comparison_peak[1] > peak[1]:
                    maximum = False
            if maximum:
                maxima.append(peak)
                appendable = False
        #Also a parameter
        elif peak[1] <= 0.01:
            appendable = True
    return(maxima)

#Open the file with the unbound spectra and read the contents
input_file = open("o_1_1.xy", "r")
input_lines = input_file.readlines()
input_file.close()
input_lines = remove_breaks(input_lines)
peaks_dict = {}
#Process each line of the text file and store the data in the lists
for line in input_lines:
    current_mz, current_intensity = line.split()
    current_mz = float(current_mz)
    current_intensity = float(current_intensity)
    if current_mz in peaks_dict.keys():
        peaks_dict[current_mz] += current_intensity
    else:
        peaks_dict[current_mz] = current_intensity
peaks_list = []
max_i = max(peaks_dict.values())
for mz, intensity in peaks_dict.items():
    peaks_list.append([mz, intensity / max_i])

local_maxima = get_maxima(peaks_list)

ub_mass = 8565
#Calculated from the weight of pt + weight of the C6N2H14 ligand
ox_mass = 195 + 114.18
ub_peaks = []
ubo_peaks = []
for charge in range(1,20):
    #moc = mass over charge
    ub_moc = ub_mass / charge
    ubo_moc = (ub_mass + ox_mass) / charge
    ub_peak = [0, 0]
    for peak in local_maxima:
        if peak[0] >= ub_moc - 1 and peak[0] <= ub_moc + 1:
            if ub_peak[1] < peak[1]:
                ub_peak = peak
    if ub_peak != [0, 0]:
        ub_peaks.append([ub_peak, charge])
        ubo_peak = [0, 0]
        for peak in local_maxima:
            if peak[0] >= ubo_moc - 1 and peak[0] <= ubo_moc + 1:
                if ubo_peak[1] < peak[1]:
                    ubo_peak = peak
        if ubo_peak != [0,0]:
            ubo_peaks.append([ubo_peak, charge])

colours = ['g','y','m','b','k','r','c','#f05e23']

import matplotlib.pyplot as plt
for i in range(0,len(ub_peaks),1):
    item = ub_peaks[i][0]
    plt.vlines(item[0],0,item[1],colors=colours[i],label='Ub (charge = ' + str(ub_peaks[i][1]) + ')')
    item = ubo_peaks[i][0]
    plt.vlines(item[0],0,item[1],colors=colours[i],linestyles='dashed',label='Ub + O')
plt.legend()
plt.show()