import matplotlib.pyplot as plt
x = list(binding_effect.keys())
y = list(binding_effect.values())
for i in range(0,len(x),100):
    plt.vlines(x[i],0,y[i])
    print(i)
plt.show()

#import matplotlib.pyplot as plt
#for i in range(0,len(filtered_difference),1):
#    item = filtered_difference[i]
#    plt.vlines(item[0],0,item[1])
#    print(i)
#plt.show()