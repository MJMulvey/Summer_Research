#import matplotlib.pyplot as plt
#x = list(binding_effect.keys())
#y = list(binding_effect.values())
#for i in range(0,len(x),100):
#    plt.vlines(x[i],0,y[i])
#    print(i)
#plt.show()


import numpy as np
from scipy.spatial.distance import euclidean

from fastdtw import fastdtw


x = np.array([[1,1], [2,2], [3,3], [4,4], [5,5]])
y = np.array([[1,1], [2,2], [3,3], [4,4], [6,6]])
distance, path = fastdtw(x, y, dist=euclidean)
print(distance)