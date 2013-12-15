import pylab as plt

datas = plt.loadtxt("en.dat")
plt.hist(datas, bins=30 , normed=1)
plt.show()
