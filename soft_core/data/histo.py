import pylab as plt

datas = plt.loadtxt("boltzmann.dat")
plt.hist(datas, bins=30)
plt.show()
