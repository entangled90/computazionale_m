import pylab as plt

datas = plt.loadtxt("tempi.dat")
plt.hist(datas, bins=1000)
plt.show()
