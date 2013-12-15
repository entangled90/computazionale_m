import pylab as plt

datas = plt.loadtxt("v2.dat")
plt.hist(datas, bins = 100,normed = 1)
plt.show()
