#Returns a dictionary of peaks obtained from taking the peaks of the original spectrum away from the peaks of the result spectrum
def get_difference(bound, unbound):
    difference_peaks = {}
    #TODO Resolution - does get_peaks provide enough? Too much?
    #TODO How should I account for small variations / noise? e.g. Two peaks in result and original with the same intensity but a tiny difference in mz (which may be the same molecule in reality)
    for mz, intensity in bound.items():
        if mz in difference_peaks:
            difference_peaks[mz] += intensity
        else:
            difference_peaks[mz] = intensity
    for mz, intensity in unbound.items():
        if mz in difference_peaks:
            difference_peaks[mz] -= intensity
        else:
            difference_peaks[mz] = -intensity
    return(difference_peaks)