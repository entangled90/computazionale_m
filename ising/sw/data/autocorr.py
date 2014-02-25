import numpy as np
import matplotlib.pyplot as plt

# generate some data
y = np.loadtxt('en_temp.dat')
# y = np.random.uniform(size=300)
result = np.zeros(100)
mean = y.mean()
for j in range(100):	
	for i in range(len(y)-j):
		result[j] += y[i]*y[i+j] / (len(y)-j)
result = result - mean*mean
result /= y.var() 
plt.plot(result)
plt.show()