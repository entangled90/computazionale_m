import pylab as plt

datas = plt.loadtxt("mean_path.dat")
plt.hist(datas, bins=50 , normed=1)
plt.show()
