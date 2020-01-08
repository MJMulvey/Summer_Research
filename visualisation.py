import matplotlib.pyplot as plt
x = list(binding_effect.keys())
y = list(binding_effect.values())
for i in range(0,len(x),1):
    plt.vlines(x[i],0,y[i])
    print(i)
plt.show()