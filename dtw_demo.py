from scipy.spatial.distance import euclidean
from fastdtw import fastdtw

s1 = [[1,1], [2,1], [3,2], [4,3], [5,2], [6,0]]
s2 = [[1,0], [2,1], [3,1], [4,2], [5,3], [6,2], [7,3], [8,3]]
s3 = [[1,0], [2,1], [3,1], [4,2], [5,3], [6,2], [7,3], [8,3]]

distance, path = fastdtw(s1, s2)
norm_distance = distance / len(s1)
print("Dist", distance)
]print(path)