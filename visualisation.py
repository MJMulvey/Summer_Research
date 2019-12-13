import matplotlib.pyplot as plt
x = list(binding_effect.keys())
y = list(binding_effect.values())
for i in range(0,len(x),100):
    print(len(x))
    plt.hlines(y[i],0,x[i])
    print(i)
plt.show()