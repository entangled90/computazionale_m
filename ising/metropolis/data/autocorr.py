import numpy as np
import matplotlib.pyplot as plt
import sys 
# generate some data

filename = sys.argv[1]

y = np.loadtxt(filename)
# y = np.random.uniform(size=300)
corr_max = 50
result = np.zeros(corr_max)
mean = y.mean()
for j in range(corr_max):	
	for i in range(len(y)-j):
		result[j] += y[i]*y[i+j] / (len(y)-j)
result = result - mean*mean
result /= y.var() 
plt.plot(result)
plt.show()